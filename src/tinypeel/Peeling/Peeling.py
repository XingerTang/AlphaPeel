from numba import jit, float32
import numpy as np


# Defining variables for peel up and peel down.
# Ideally these would be characters, but numba does not support characters.
PEEL_UP = 0
PEEL_DOWN = 1


@jit(
    nopython=True,
    nogil=True,
    locals={"e": float32, "e4": float32, "e16": float32, "e1e": float32},
    fastmath=True,
)
def peel(family, operation, peelingInfo, singleLocusMode):
    """
    Main peeling function.
    """
    isSexChrom = peelingInfo.isSexChrom

    e = 0.000001
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
    jointParents = jointParents * e1e + e / 16
    # There are 4 x 4 values for each locus in jointparents.

    # There are 4 values for each locus. Normalization is done here so that jointParents is as accurate as possible.
    # We need the posterior terms here for the peeling up step later.
    probSire = probSire * e1e + e4
    probDam = probDam * e1e + e4

    # Now construct the parental genotypes based on within-family information.

    for index in range(nOffspring):
        child = family.offspring[index]

        # The child's estimate is the combination of the posterior term and penetrance term for that child.
        # We are estimating the parent's genotypes so the anterior term is ignored to avoid double counting.
        childValues = posterior[child] * penetrance[child]
        childValues /= summing(childValues)
        childValues = e1e * childValues + e4

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
        # This is done to take away the setimate for an individual child from their parent's posterior term.
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
        sirePosterior = sirePosterior * e1e + e4
        peelingInfo.posteriorSire_new[fam] = sirePosterior

        damPosterior = combineAndReduceAxis0(allToParents, probSire)
        damPosterior /= summing(damPosterior)
        damPosterior = damPosterior * e1e + e4
        peelingInfo.posteriorDam_new[fam] = damPosterior

    if (not singleLocusMode) and (operation == PEEL_DOWN):
        # Estimate the segregation probabilities for each child.

        for i in range(nOffspring):
            # Child values is the same as in the posterior estimation step above.
            child = family.offspring[i]
            childValues = posterior[child] * penetrance[child]
            childValues /= summing(childValues)
            childValues = e1e * childValues + e4

            if isSexChrom and peelingInfo.sex[child] == 0:  # 0 for male, 1 for female.
                segregationTensor = peelingInfo.segregationTensorXY
            if isSexChrom and peelingInfo.sex[child] == 1:  # 0 for male, 1 for female.
                segregationTensor = peelingInfo.segregationTensorXX

            # Option 1: Estimate without normalizing.
            # estimateSegregation(segregationTensor, parentsMinusChild[i,:,:,:], childValues, pointSeg[child,:,:])
            # Option 2: Estimate with normalizing. I think this is what we want.
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
            )


# The following are a large number of "helper" jit functions that replace the einstien sums in the original scripts.


@jit(nopython=True, nogil=True, fastmath=True)
def getJointParents(probSire, probDam):
    """
    Create the joint parental genotypes based on the probabilities for each parent.
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
    """
    Create the child-specific segregation tensor using the child's currrent segregation estimate.
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
    """
    Estimate the parental genotypes based on the child's genotypes and their segregation tensor.
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
    """
    Project the parent genotypes down onto the child genotypes
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
    """
    Estimate with normalizing.
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
    """
    Matrix is 4 x 4 x nLoci:
    Output is to take the exponential of the matrix and normalize each locus.
    We need to make sure that there are not any overflow values.
    """
    nLoci = mat.shape[2]
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
    """
    Matrix is 4 x 4 x nLoci
    Output is to take the exponential of the matrix and normalize each locus.
    We need to make sure that there are not any overflow values.
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
    """
    This is the forward backward algorithm.
    Segregation estimate state ordering: pp, pm, mp, mm
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
        e = transmission[locus - 1]
        e2 = e**2
        e1e = e * (1.0 - e)
        e2i = (1.0 - e) ** 2
        for state in range(4):
            tmp[state] = prev[state] * pointSeg[state, locus - 1]

        sum_state = 0.0
        for state in range(4):
            sum_state += tmp[state]
        for state in range(4):
            tmp[state] = tmp[state] / sum_state

        # !                  fm  fm  fm  fm
        # !segregationOrder: pp, pm, mp, mm
        new[0] = e2 * tmp[3] + e1e * (tmp[1] + tmp[2]) + e2i * tmp[0]
        new[1] = e2 * tmp[2] + e1e * (tmp[0] + tmp[3]) + e2i * tmp[1]
        new[2] = e2 * tmp[1] + e1e * (tmp[0] + tmp[3]) + e2i * tmp[2]
        new[3] = e2 * tmp[0] + e1e * (tmp[1] + tmp[2]) + e2i * tmp[3]

        for state in range(4):
            seg[state, locus] *= new[state]
        prev = new

    prev = np.full((4), 0.25, dtype=np.float32)
    for locus in range(
        nLoci - 2, -1, -1
    ):  # zero indexed then minus one since we skip the boundary.
        e = transmission[locus]
        e2 = e**2
        e1e = e * (1.0 - e)
        e2i = (1.0 - e) ** 2

        for state in range(4):
            tmp[state] = prev[state] * pointSeg[state, locus + 1]

        sum_state = 0.0
        for state in range(4):
            sum_state += tmp[state]
        for state in range(4):
            tmp[state] = tmp[state] / sum_state

        new[0] = e2 * tmp[3] + e1e * (tmp[1] + tmp[2]) + e2i * tmp[0]
        new[1] = e2 * tmp[2] + e1e * (tmp[0] + tmp[3]) + e2i * tmp[1]
        new[2] = e2 * tmp[1] + e1e * (tmp[0] + tmp[3]) + e2i * tmp[2]
        new[3] = e2 * tmp[0] + e1e * (tmp[1] + tmp[2]) + e2i * tmp[3]

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
    """
    Summing over axis 0,
    equivalent to np.sum(vec, axis=0)
    """
    total = np.zeros(vec.shape[1:], dtype=np.float32)
    for i in range(vec.shape[0]):
        total += vec[i]
    return total


@jit(nopython=True, nogil=True, fastmath=True)
def summing_twice(vec):
    """
    Summing over axis 0 and 1,
    equivalent to np.sum(vec, axis=(0, 1))
    """
    total = np.zeros(vec.shape[2:], dtype=np.float32)
    for i in range(vec.shape[0]):
        for j in range(vec.shape[1]):
            total += vec[i, j]
    return total
