from numba import jit, float32
import numpy as np

PEEL_UP = 0
PEEL_DOWN = 1


@jit(
    nopython=True,
    nogil=True,
    locals={"e": float32, "e4": float32, "e1e": float32},
    fastmath=True,
)
def peel(family, operation, peelingInfo, singleLocusMode):
    """This is the main peeling function.

    :param family: The family object that the peeling is performed on.
    :type family: class:`tinyhouse.Pedigree.Family`
    :param operation: A flag to indicate the direction of the peeling process.
        `0` if peeling up, and `1` if peeling down.
    :type operation: int
    :param peelingInfo: Peeling information container.
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :param singleLocusMode: A flag to indicate the mode of peeling.
        `False` if using multi-locus peeling, and `True` if using single-locus peeling.
    :type singleLocusMode: bool
    """
    isSexChrom = peelingInfo.isSexChrom

    # An error term `e` and associated values that are used frequently
    e = 0.00106  # First definition in peel function?
    e1e = 1.0 - e
    e4 = e / 4

    # Setup local variables from the peeling information container.
    anterior = peelingInfo.anterior
    penetrance = peelingInfo.penetrance
    posterior = peelingInfo.posterior
    segregation = peelingInfo.segregation

    pointSeg = peelingInfo.pointSeg
    segregationTensor = peelingInfo.segregationTensor
    segregationTensor_norm = peelingInfo.segregationTensor_norm

    nLoci = peelingInfo.nLoci
    nOffspring = len(family.offspring)
    sire = family.sire
    dam = family.dam
    fam = family.idn

    # Creating variables here:
    # childToParents: The projection of each child onto the parent genotypes.
    # childSegs: The segregation estimates for a particular child (These are re-used? so need to be stored)
    # allToParents: The projection of each child onto the parental genotypes.
    # parentsMinustChild: The estimate of the parent's genotypes minus the contribution from a specific individual.

    childToParents = np.zeros((nOffspring, 4, 4, nLoci), dtype=np.float32)
    childSegTensor = np.zeros((nOffspring, 4, 4, 4, nLoci), dtype=np.float32)
    allToParents = np.ones((4, 4, nLoci), dtype=np.float32)
    parentsMinusChild = np.ones((nOffspring, 4, 4, nLoci), dtype=np.float32)

    # Some local variables. currentSeg is the segregation estimate of a child (but may be modified).
    currentSeg = np.ones((4, nLoci), dtype=np.float32)

    # Construct the joint parent genotypes based on the parent's anterior, penetrance, and posterior terms minus this family.
    probSire = (
        anterior[sire] * penetrance[sire] * peelingInfo.posteriorSire_minusFam[fam]
    )
    probDam = anterior[dam] * penetrance[dam] * peelingInfo.posteriorDam_minusFam[fam]

    probSire /= summing(probSire)
    probDam /= summing(probDam)

    # Create the joint parental genotypes based on the probabilities for each parent.
    jointParents = getJointParents(probSire, probDam)
    jointParents /= summing_twice(jointParents)
    jointParents = jointParents * e1e + e / 16  # Unclear error for Jointparents
    # There are 4 x 4 values for each locus in jointparents.

    # There are 4 values for each locus. Normalization is done here so that jointParents is as accurate as possible.
    # We need the posterior terms here for the peeling up step later.
    probSire = probSire * e1e + e4  # Unclear error for probSire
    probDam = probDam * e1e + e4  # Unclear error for probDam

    # Now construct the parental genotypes based on within-family information.

    for index in range(nOffspring):
        child = family.offspring[index]

        # The child's estimate is the combination of the posterior term and penetrance term for that child.
        # We are estimating the parent's genotypes so the anterior term is ignored to avoid double counting.
        childValues = posterior[child] * penetrance[child]
        childValues /= summing(childValues)
        childValues = e1e * childValues + e4  # Unclear error for childValues

        # METHOD 1: Just use the current segregation of the child.
        currentSeg[:, :] = segregation[child, :, :]
        currentSeg /= summing(currentSeg)

        if isSexChrom and peelingInfo.sex[child] == 0:  # 0 for male, 1 for female.
            segregationTensor = peelingInfo.segregationTensorXY
        if isSexChrom and peelingInfo.sex[child] == 1:  # 0 for male, 1 for female.
            segregationTensor = peelingInfo.segregationTensorXX

        # Create the child-specific segregation tensor using the child's currrent segregation estimate.
        createChildSegs(segregationTensor, currentSeg, childSegTensor[index])

        # Estimate the parental genotypes based on the child's genotypes and their segregation tensor.
        projectChildGenotypes(
            childSegTensor[index],
            childValues,
            childToParents[index],
        )

    # Estimate the parents genotype and the child-specific posterior terms using a slightly smarter log scale.
    allToParents[:, :, :] = 0.0
    for i in range(nOffspring):
        parentsMinusChild[i] = np.log(jointParents)
        log_childToParents = np.log(childToParents[i])
        allToParents += log_childToParents
        parentsMinusChild[i] -= log_childToParents
        # This is done to take away the estimate for an individual child from their parent's posterior term.
    parentsMinusChild += allToParents

    # Move from a log-scale to a non-log scale and re-normalize.
    allToParents = expNorm2D(allToParents)
    for i in range(nOffspring):
        parentsMinusChild[i] = expNorm2D(parentsMinusChild[i])

    if operation == PEEL_DOWN:
        for i in range(nOffspring):
            child = family.offspring[i]

            # Project the parent genotypes down onto the child genotypes
            projectParentGenotypes(
                childSegTensor[i],
                parentsMinusChild[i],
                anterior[child],
            )
            anterior[child] /= summing(anterior[child])

    if operation == PEEL_UP:
        # Take the allToParents estimate and combine to
        # estimate the sire and dam's posterior estimates (for this family)

        sirePosterior = combineAndReduceAxis1(allToParents, probDam)
        sirePosterior /= summing(sirePosterior)
        sirePosterior = sirePosterior * e1e + e4  # Unclear error for sirePosterior
        peelingInfo.posteriorSire_new[fam] = sirePosterior

        damPosterior = combineAndReduceAxis0(allToParents, probSire)
        damPosterior /= summing(damPosterior)
        damPosterior = damPosterior * e1e + e4  # Unclear error for damPosterior
        peelingInfo.posteriorDam_new[fam] = damPosterior

    if (not singleLocusMode) and (operation == PEEL_DOWN):
        # Estimate the segregation probabilities for each child.

        for i in range(nOffspring):
            # Child values is the same as in the posterior estimation step above.
            child = family.offspring[i]
            childValues = posterior[child] * penetrance[child]
            childValues /= summing(childValues)
            childValues = (
                e1e * childValues + e4
            )  # Unclear error for childValues during peel down

            if isSexChrom and peelingInfo.sex[child] == 0:  # 0 for male, 1 for female.
                segregationTensor = peelingInfo.segregationTensorXY
            if isSexChrom and peelingInfo.sex[child] == 1:  # 0 for male, 1 for female.
                segregationTensor = peelingInfo.segregationTensorXX

            # Estimate with normalizing.
            estimateSegregationWithNorm(
                segregationTensor,
                segregationTensor_norm,
                parentsMinusChild[i],
                childValues,
                pointSeg[child],
            )

            segregation[child] = (
                e1e * collapsePointSeg(pointSeg[child], peelingInfo.transmissionRate)
                + e4
            )  # Unclear error for segregation during peel down


# The following are a large number of "helper" jit functions that replace the einstien sums in the original scripts.


@jit(nopython=True, nogil=True, fastmath=True)
def getJointParents(probSire, probDam):
    """Create the joint parental genotypes based on the probabilities for each parent.

    :param probSire: the probability of each genotype of each locus of the sire
        with information of the current child from previous peeling cycle (P(p))
    :type probSire: 2D numpy array of float32 with size 4 x nLoci
    :param probDam: the probability of each genotype of each locus of the dam
        with information of the current child from previous peeling cycle (P(m))
    :type probDam: 2D numpy array of float32 with size 4 x nLoci
    :return: the probabilities of all the combinations of the sire's genotype and the dam's genotype
        with information from previous peeling cycle (P(p, m))
    :rtype: 3D numpy array of float32 with size 4 x 4 x nLoci
    """
    nLoci = probSire.shape[1]
    output = np.zeros(shape=(4, 4, nLoci), dtype=np.float32)
    for p in range(4):
        for m in range(4):
            for locus in range(nLoci):
                output[p, m, locus] = probSire[p, locus] * probDam[m, locus]
    return output


@jit(nopython=True, nogil=True, fastmath=True)
def createChildSegs(segregationTensor, currentSeg, output):
    """Create the child-specific segregation tensor using the child's current segregation estimate.

    :param segregationTensor: the probability of each combination of the sire's genotype, the dam's genotype,
        child's genotype and segregation without any other information (P(p, m, allele, seg))
    :type segregationTensor: 4D numpy array of float32 with size 4 x 4 x 4 x 4
    :param currentSeg: The probability of each segregation of each locus of the child
        with information from previous peeling cycle (P(seg))
    :type currentSeg: 2D numpy array of float32 with size 4 x nLoci
    :param output: the probability of each combination of the sire's genotype, the dam's genotype and
        the child's genotype of each locus with information from previous peeling cycle (P(p, m, allele))
    :type output: 4D numpy array of float32 with size 4 x 4 x 4 x nLoci
    """
    nLoci = currentSeg.shape[1]
    output[:, :, :, :] = 0.0
    for p in range(4):
        for m in range(4):
            for allele in range(4):
                for seg in range(4):
                    for locus in range(nLoci):
                        output[p, m, allele, locus] += (
                            segregationTensor[p, m, allele, seg]
                            * currentSeg[seg, locus]
                        )


@jit(nopython=True, nogil=True, fastmath=True)
def projectChildGenotypes(childSegs, childValues, output):
    """Estimate the parental genotypes based on the child's genotypes and their segregation tensor.

    :param childSegs: the probability of each combination of the sire's genotype, the dam's genotype and
        the child's genotype of each locus with information from previous peeling cycle
        (P(p, m, allele))
    :type childSegs: 4D numpy array of float32 with size 4 x 4 x 4 x nLoci
    :param childValues: the probability of each genotype of each locus of the current child
        given the information of itself and its later generations from previous peeling cycle
        (P(allele))
    :type childValues: 2D numpy array of float32 with size 4 x nLoci
    :param output: the probability of each combination of the sire's genotype, the dam's genotype
        given the information of later and current generations from previous peeling cycle
        (P(p, m))
    :type output: 3D numpy array of float32 with size 4 x 4 x nLoci
    """
    nLoci = childSegs.shape[3]
    output[:, :, :] = 0.0
    for p in range(4):
        for m in range(4):
            for allele in range(4):
                for locus in range(nLoci):
                    output[p, m, locus] += (
                        childSegs[p, m, allele, locus] * childValues[allele, locus]
                    )


@jit(nopython=True, nogil=True, fastmath=True)
def projectParentGenotypes(childSegs, parentValues, output):
    """Project the parent genotypes down onto the child genotypes.

    :param childSegs: the probability of each combination of the sire's genotype, the dam's genotype and
        the child's genotype of each locus with information from previous peeling cycle
        (P(p, m, allele))
    :type childSegs: 4D numpy array of float32 with size 4 x 4 x 4 x nLoci
    :param parentValues: the probability of each combination of sire's genotype and the dams's genotype
        of each locus given the information of later and current generations from previous peeling cycle
        without the information of the current child (P(p, m))
    :type parentValues: 3D numpy array of float32 with size 4 x 4 x nLoci
    :param output: the probability of child's genotype of each locus given the later and current generations
        from previous peeling cycle without the information of the current child (P(allele))
    :type output: 2D numpy array of float32 with size 4 x nLoci
    """
    nLoci = childSegs.shape[3]
    output[:, :] = 0.0
    for p in range(4):
        for m in range(4):
            for allele in range(4):
                for locus in range(nLoci):
                    output[allele, locus] += (
                        childSegs[p, m, allele, locus] * parentValues[p, m, locus]
                    )


@jit(nopython=True, nogil=True, fastmath=True)
def estimateSegregationWithNorm(
    segregationTensor, segregationTensor_norm, parentValues, childValues, output
):
    """Estimate with normalizing.

    :param segregationTensor: the probability of each combination of the sire's genotype, the dam's genotype and
        the child's genotype and segregation without any other information (P(p, m, allele, seg))
    :type segregationTensor: 4D numpy array of float32 with size 4 x 4 x 4 x 4
    :param segregationTensor_norm: the mean probability of each combination of the sire's genotype, the dam's genotype
        and the child's genotype across the child's segregation without any other information
        (1 / 4 x (P(p, m, allele)))
    :type segregationTensor_norm: 3D numpy array of float32 with size 4 x 4 x 4
    :param parentValues: the probability of each combination of sire's genotype and the dams's genotype
        of each locus given the information of later and current generations from previous peeling cycle
        without the information of the current child (P(p, m))
    :type parentValues: 3D numpy array of float32 with size 4 x 4 x nLoci
    :param childValues: the probability of each genotype of each locus of the current child
        given the information of itself and its later generations from previous peeling cycle
        (P(allele))
    :type childValues: 2D numpy array of float32 with size 4 x nLoci
    :param output: the probability of each segregation states of each locus of the current child
        given the information of later and current generations from previous peeling cycle
        (P(seg))
    :type output: 2D numpy array of float32 with size 4 x nLoci
    """
    nLoci = childValues.shape[1]
    output[:, :] = 0.0
    for p in range(4):
        for m in range(4):
            for allele in range(4):
                for seg in range(4):
                    for locus in range(nLoci):
                        # Check if norm is 0. Otherwise use norm to normalize.
                        if segregationTensor_norm[p, m, allele] != 0:
                            output[seg, locus] += (
                                segregationTensor[p, m, allele, seg]
                                * parentValues[p, m, locus]
                                * childValues[allele, locus]
                                / segregationTensor_norm[p, m, allele]
                            )


@jit(nopython=True, nogil=True, fastmath=True)
def combineAndReduceAxis1(jointEstimate, parentEstimate):
    """Summing over axis 1 of jointEstimate with weights given by parentEstimate

    :param jointEstimate: the probability of each combination of sire's genotype and the dams's genotype
        of each locus given the information of later and current generations from previous peeling cycle
        (P(p, m))
    :type jointEstimate: 3D numpy array of float32 with size 4 x 4 x nLoci
    :param parentEstimate: the probability of each genotype of each locus of the dam
        with information of the current child from previous peeling cycle (P(m))
    :type parentEstimate: 2D numpy array of float32 with size 4 x nLoci
    :return: the probability of each genotype of each locus of the sire
        given the information of later and current generations from previous peeling cycle (P(p))
    :rtype: 2D numpy array of float32 with size 4 x nLoci
    """
    nLoci = parentEstimate.shape[1]
    output = np.zeros((4, nLoci), dtype=np.float32)
    for p in range(4):
        for m in range(4):
            for locus in range(nLoci):
                output[p, locus] += (
                    jointEstimate[p, m, locus] * parentEstimate[m, locus]
                )
    return output


@jit(nopython=True, nogil=True, fastmath=True)
def combineAndReduceAxis0(jointEstimate, parentEstimate):
    """Summing over axis 0 of jointEstimate with weights given by parentEstimate

    :param jointEstimate: the probability of each combination of sire's genotype and the dams's genotype
        of each locus given the information of later and current generations from previous peeling cycle
        (P(p, m))
    :type jointEstimate: 3D numpy array of float32 with size 4 x 4 x nLoci
    :param parentEstimate: the probability of each genotype of each locus of the sire
        with information of the current child from previous peeling cycle (P(p))
    :type parentEstimate: 2D numpy array of float32 with size 4 x nLoci
    :return: the probability of each genotype of each locus of the dam
        given the information of later and current generations from previous peeling cycle (P(m))
    :rtype: 2D numpy array of float32 with size 4 x nLoci
    """
    nLoci = parentEstimate.shape[1]
    output = np.zeros((4, nLoci), dtype=np.float32)
    for p in range(4):
        for m in range(4):
            for locus in range(nLoci):
                output[m, locus] += (
                    jointEstimate[p, m, locus] * parentEstimate[p, locus]
                )
    return output


@jit(nopython=True, nogil=True, fastmath=True)
def expNorm2D(mat):
    """Output is to take the exponential of the matrix and normalize each locus.

    :param mat: a 3D matrix with last axis represents the locus
    :type mat: 3D numpy array of float32 with size 4 x 4 x nLoci
    :return: the normalized exponential of the `mat`
    :rtype: 3D numpy array of float32 with size 4 x 4 x nLoci
    """
    nLoci = mat.shape[2]
    # We need to make sure that there are not any overflow values.
    for locus in range(nLoci):
        maxVal = mat[0, 0, locus]
        for p in range(4):
            for m in range(4):
                if mat[p, m, locus] > maxVal:
                    maxVal = mat[p, m, locus]
        for p in range(4):
            for m in range(4):
                mat[p, m, locus] -= maxVal
    # Normalize.
    tmp = np.exp(mat)
    for locus in range(nLoci):
        total = 0
        for p in range(4):
            for m in range(4):
                total += tmp[p, m, locus]
        for p in range(4):
            for m in range(4):
                tmp[p, m, locus] /= total
    return tmp


@jit(nopython=True, nogil=True, fastmath=True)
def expNorm1D(mat):
    """Output is to take the exponential of the matrix and normalize each locus.

    :param mat: a 2D matrix with last axis represents the locus
    :type mat: 2D numpy array of float32 with size 4 x nLoci
    :return: the normalized exponential of the `mat`
    :rtype: 2D numpy array of float32 with size 4 x nLoci
    """

    nLoci = mat.shape[1]
    for locus in range(nLoci):
        maxVal = mat[0, locus]
        for p in range(4):
            if mat[p, locus] > maxVal:
                maxVal = mat[p, locus]
        for p in range(4):
            mat[p, locus] -= maxVal
    tmp = np.exp(mat)
    for locus in range(nLoci):
        total = 0
        for p in range(4):
            total += tmp[p, locus]
        for p in range(4):
            tmp[p, locus] /= total
    return tmp


@jit(
    nopython=True,
    nogil=True,
    locals={"e": float32, "e2": float32, "e1e": float32, "e2i": float32},
    fastmath=True,
)
def collapsePointSeg(pointSeg, transmission):
    """Using Baum-Welch algorithm to calculate the segregation probabilities.

    :param pointSeg: the probability of each segregation states of each locus of the current child
        given the information of later and current generations from previous peeling cycle
        (P(seg))
        the segregation state ordering: pp, pm, mp, mm
    :type pointSeg: 2D numpy array of float32 with size 4 x nLoci
    :param transmission: transmission function based on the distance
        (possibly recombination rate in future)
    :type transmission: 1D numpy array of float32 with size (nLoci - 1)
    :return: the probability of each segregation states of each locus of the current child
        after the implemtation of Baum-Welch algorithm
    """
    nLoci = pointSeg.shape[1]

    seg = np.full(pointSeg.shape, 0.25, dtype=np.float32)
    for locus in range(nLoci):
        for state in range(4):
            seg[state, locus] = pointSeg[state, locus]

    tmp = np.zeros((4), dtype=np.float32)
    new = np.zeros((4), dtype=np.float32)

    prev = np.full((4), 0.25, dtype=np.float32)
    for locus in range(1, nLoci):
        r = transmission[locus - 1]
        r2 = r**2
        r1r = r * (1.0 - r)
        r2i = (1.0 - r) ** 2
        for state in range(4):
            tmp[state] = prev[state] * pointSeg[state, locus - 1]

        sum_state = 0.0
        for state in range(4):
            sum_state += tmp[state]
        for state in range(4):
            tmp[state] = tmp[state] / sum_state

        # !                  fm  fm  fm  fm
        # !segregationOrder: pp, pm, mp, mm
        new[0] = r2 * tmp[3] + r1r * (tmp[1] + tmp[2]) + r2i * tmp[0]
        new[1] = r2 * tmp[2] + r1r * (tmp[0] + tmp[3]) + r2i * tmp[1]
        new[2] = r2 * tmp[1] + r1r * (tmp[0] + tmp[3]) + r2i * tmp[2]
        new[3] = r2 * tmp[0] + r1r * (tmp[1] + tmp[2]) + r2i * tmp[3]

        for state in range(4):
            seg[state, locus] *= new[state]
        prev = new

    prev = np.full((4), 0.25, dtype=np.float32)
    for locus in range(
        nLoci - 2, -1, -1
    ):  # zero indexed then minus one since we skip the boundary.
        r = transmission[locus]
        r2 = r**2
        r1r = r * (1.0 - r)
        r2i = (1.0 - r) ** 2

        for state in range(4):
            tmp[state] = prev[state] * pointSeg[state, locus + 1]

        sum_state = 0.0
        for state in range(4):
            sum_state += tmp[state]
        for state in range(4):
            tmp[state] = tmp[state] / sum_state

        new[0] = r2 * tmp[3] + r1r * (tmp[1] + tmp[2]) + r2i * tmp[0]
        new[1] = r2 * tmp[2] + r1r * (tmp[0] + tmp[3]) + r2i * tmp[1]
        new[2] = r2 * tmp[1] + r1r * (tmp[0] + tmp[3]) + r2i * tmp[2]
        new[3] = r2 * tmp[0] + r1r * (tmp[1] + tmp[2]) + r2i * tmp[3]

        for state in range(4):
            seg[state, locus] *= new[state]
        prev = new

    for locus in range(nLoci):
        sum_state = 0.0
        for state in range(4):
            sum_state += seg[state, locus]
        for state in range(4):
            seg[state, locus] = seg[state, locus] / sum_state

    return seg


@jit(nopython=True, nogil=True, fastmath=True)
def summing(vec):
    """Summing the input `vec` over axis 0, equivalent to np.sum(vec, axis=0).

    :param vec: a vector
    :type vec: numpy array with dimension > 1
    :return: a vector
    :rtype: numpy array with one less dimension compared to `vec`
    """
    total = np.zeros(vec.shape[1:], dtype=np.float32)
    for i in range(vec.shape[0]):
        total += vec[i]
    return total


@jit(nopython=True, nogil=True, fastmath=True)
def summing_twice(vec):
    """Summing the input `vec` over axis 0 and 1, equivalent to np.sum(vec, axis=(0, 1)).

    :param vec: a vector
    :type vec: numpy array with dimension > 2
    :return: a vector
    :rtype: numpy array with one less dimension compared to `vec`
    """
    total = np.zeros(vec.shape[2:], dtype=np.float32)
    for i in range(vec.shape[0]):
        for j in range(vec.shape[1]):
            total += vec[i, j]
    return total
