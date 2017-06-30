def construct_graph(filename,symmetry=xtal.Symmetry('Cubic')):
    fzQu = fuZqu(symmetry)
    f = h5py.File(filename, 'r')['DataContainers']
    try:
        data = f['ImageDataContainer']
    except KeyError:
        try:
            data = f['SyntheticVolumeContainer']
        except KeyError:
            raise KeyError('Dream3D File must contain ImageDataContainer or SyntheticVolumeContainer')

    interior_grains = data['CellFeatureData']['NumCells'][:] * (1-data['CellFeatureData']['SurfaceFeatures'])
    grainIDs = np.argsort(voxels_grains[:,0])
    
    #considering only sigma3 boundaries as twin boundaries
    misorientation_list = np.array(data['CellFeatureData']['MisorientationList'])
    twin_gb_list = np.array(list(map(is_sigma_3,misorientation_list)))

    neighbor_list = np.array(data['CellFeatureData']['NeighborList'])
    num_neighbors = np.array(data['CellFeatureData']['NumNeighbors'])
    node_indices = np.cumsum(num_neighbors+1)

    voxels_interior = f['DataContainers/ImageDataContainer/CellFeatureData/NumCells'][:] * (1-f['DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures'][:])
    clusterId = np.argsort(voxels_interior[:,0])
    internalClusters = [i for i in clusterId if voxels_interior[i] > 0]

    clusterCount = 0
    with open('cluster_stats2.csv', 'w') as outfile:
        outfile.write('clusterId, numGrains, numVariants, twinArea, totalArea, clusterVolume, maxVolume, minDist, maxDist, numGenerations, numLines, numVarLines\n')
        for clusterId in internalClusters:
        # # clusterId = clusterId[-2]
        # clusterId = clusterId[-50]
        # clusterId = 769
        # print(clusterId)

            # get grains belonging to cluster
            clusterGrains = np.nonzero(f['DataContainers/ImageDataContainer/CellFeatureData/ClusterIds'][:]==clusterId)[0]
            if 1 is len(clusterGrains):
                continue
            clusterCount += 1

            # get orientations of grains in cluster and convert from xyzw to wxyz
            quats = np.roll(f['DataContainers/ImageDataContainer/CellFeatureData/AvgQuats'][:][clusterGrains], 1, axis = 1)

            # get sizes of grains in cluster
            diams = f['DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters'][:][clusterGrains]

            # get neighborss of grains in cluster
            neighbors = f['DataContainers/ImageDataContainer/CellFeatureData/NeighborList']
            neighborCount = f['DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors'][:,0]
            neighborInds = np.copy(neighborCount)
            for i in range(1, len(neighborInds)):
                neighborInds[i] += neighborInds[i-1]
            neighborInds -= neighborCount

            # get shared surface areas
            sharedAreas = f['DataContainers/ImageDataContainer/CellFeatureData/SharedSurfaceAreaList']

            # build list of neighbors within cluster
            lines = []
            areas = []
            for i in range(clusterGrains.shape[0]):
                grain = clusterGrains[i]
                offset = neighborInds[grain]
                count = neighborCount[grain]
                grainNeighbors = neighbors[offset:offset+count]
                for n in range(len(grainNeighbors)):
                    ind = np.argwhere(clusterGrains == grainNeighbors[n])[:,0]
                    if ind.shape[0] > 0:
                        lines.append((i,ind[0]))
                        areas.append(sharedAreas[offset+n])

            # segment cluster into variants
            variants = np.zeros_like(clusterGrains)
            variantId = 1
            for i in range(clusterGrains.shape[0]):
                for j in range(i):
                    miso = Symmetry.Cubic.disoQuat(quats[i], quats[j])
                    if 2*np.math.acos(min(miso[0],1))*180/np.pi <= 5: # within 5 degrees
                        variants[i] = variants[j]
                        break
                if 0 == variants[i]:
                    variants[i] = variantId
                    variantId += 1
            variantId -= 1
            variants -= 1

            # move quats so that they are nearby eachother
            fzQuats = [Symmetry.Cubic.fzQu(Quaternion(q)) for q in quats]
            baseQuats = [None] * variantId
            for i in range(len(variants)):
                for j in range(variantId):
                    if baseQuats[variants[j]] is None:
                        baseQuats[variants[j]] = Quaternion(fzQuats[i])
                    else:
                        fzQuats[i] = Symmetry.Cubic.nearbyQuat(baseQuats[variants[j]], fzQuats[i])

            # compute volume weighted average orientation of each variant
            grainCus = np.array([rotations.qu2cu(Symmetry.Cubic.fzQu(q.wxyz)) for q in fzQuats])

            avgCu = np.zeros((variantId, 3))
            vol = np.zeros((variantId,))
            for i in range(clusterGrains.shape[0]):
                v = variants[i]
                avgCu[v] += grainCus[i,:] * (diams[i]**3)
                vol[v] += diams[i]**3

            avgCu[:,:] /= vol[:,np.newaxis]

            avgQu = []
            for cu in avgCu:
                avgQu.append(Quaternion(rotations.cu2qu(list(cu))))

            # build connections between variants
            variantLines = []
            variantAreas = []
            lineCounts = []
            for i in range(len(lines)):
                line = lines[i]
                variantPair = sorted((variants[line[0]], variants[line[1]]))
                if variantPair[0] != variantPair[1]:
                    if variantPair in variantLines:
                        ind = variantLines.index(variantPair)
                        variantAreas[ind] += areas[i]
                        lineCounts[ind] += 1
                    else:
                        variantLines.append(variantPair)
                        variantAreas.append(0)
                        lineCounts.append(1)
            variantAreas = np.sqrt(variantAreas)

            variantMisos = []
            twinLine = []
            twinLines = []
            for line in variantLines:
                miso = Symmetry.Cubic.disoQuat(avgQu[line[0]], avgQu[line[1]])
                variantMisos.append(miso)
                nDot = sum(miso[1:]) / (np.math.sqrt(3) * np.math.sqrt(1 - miso[0]*miso[0]))
                rot_60 = abs(2*np.math.acos(miso[0])-np.pi/3) * 180 / np.pi <= 10 # within 5 degrees of 60 rotation
                axis_111 = np.math.acos(min(nDot, 1)) * 180 / np.pi <= 10 # within 5 degrees of 111 axis
                if rot_60 and axis_111:
                    twinLine.append(True)
                    twinLines.append(line)
                else:
                    twinLine.append(False)

            # find largest variant
            varMax = np.argmax(vol)
            network = UndirectedGraph(twinLines)

            # compute network metrics
            maxDist = network.maxDist()
            volRootDist = network.computeDistance(varMax)
            distRoots = network.findMinRoot()
            minDist = -1
            numGenerations = -1
            if distRoots is not None:
                distRootDist = network.computeDistance(distRoots[0])
                seperation = distRootDist[varMax]
                equivRoots = len(distRoots)

                generations = volRootDist
                # print(generations)
                numGenerations = max(generations)
                minDist = max(distRootDist)

                numGrains = len(variants) # number of grain in twin related domain
                numVariants = len(generations) # number of varaints
                # print(numGrains)
                clusterVolume = sum(vol) # total volume of twin related domain
                twinArea = 0
                for i in range(len(variantLines)):
                    if twinLine[i]:
                        twinArea +=variantAreas[i]
                totalArea = sum(variantAreas)
                # print(twinArea / totalArea)
                print(clusterId, minDist, maxDist, numGenerations)
                outfile.write('%i,%i,%i,%f,%f,%f,%f,%i,%i,%i,%f,%f\n' % (clusterId, numGrains, numVariants, twinArea, totalArea, clusterVolume, max(vol), minDist, maxDist, numGenerations, len(lines), len(variantLines)))
    print(clusterCount)