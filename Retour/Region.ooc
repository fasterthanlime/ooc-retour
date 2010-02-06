use retour
import Retour/[Retour]

calculateDistanceField: static func(chf: RCCompactHeightfield@, src, dst: UShort*, maxDist: UShort@) -> UShort* {
	w := chf width
	h := chf height
	
	// Init distance and points.
	for (i: Int in 0..(chf spanCount))
		src[i] = 0xffff
	
	// Mark boundary cells.
	for (y: Int in 0..h) {
		for (x: Int in 0..w) {
			c := chf cells[x+y*w]
			for (i: Int in (c index as Int)..((c index + c count) as Int)) {
				s := chf spans[i]
				area := chf areas[i]
				
				nc := 0
				for (dir: Int in 0..4) {
					if (rcGetCon(s, dir) != RC_NOT_CONNECTED) {
						ax := x + rcGetDirOffsetX(dir)
						ay := y + rcGetDirOffsetY(dir)
						ai := (chf cells[ax+ay*w] index as Int) + rcGetCon(s, dir)
						if (area == chf areas[ai])
							nc+=1
					}
				}
				if (nc != 4)
					src[i] = 0
			}
		}
	}
	
	// Pass 1
	for (y: Int in 0..h) {
		for (x: Int in 0..w) {
			c := chf cells[x+y*w]
			for (i: Int in (c index as Int)..(c index + c count) as Int) {
				s := chf spans[i] 
				if (rcGetCon(s, 0) != RC_NOT_CONNECTED) {
					// (-1,0)
					ax := x + rcGetDirOffsetX(0)
					ay := y + rcGetDirOffsetY(0)
					ai := (chf cells[ax+ay*w] index as Int) + rcGetCon(s, 0)
					asz := chf spans[ai]
					if (src[ai]+2 < src[i])
						src[i] = src[ai]+2
					
					// (-1,-1)
					if (rcGetCon(asz, 3) != RC_NOT_CONNECTED) {
						aax := ax + rcGetDirOffsetX(3)
						aay := ay + rcGetDirOffsetY(3)
						aai := (chf cells[aax+aay*w] index as Int) + rcGetCon(asz, 3)
						if (src[aai]+3 < src[i])
							src[i] = src[aai]+3
					}
				}
				if (rcGetCon(s, 3) != RC_NOT_CONNECTED) {
					// (0,-1)
					ax := x + rcGetDirOffsetX(3)
					ay := y + rcGetDirOffsetY(3)
					ai := (chf cells[ax+ay*w] index as Int) + rcGetCon(s, 3)
					asz := chf spans[ai]
					if (src[ai]+2 < src[i])
						src[i] = src[ai]+2
					
					// (1,-1)
					if (rcGetCon(asz, 2) != RC_NOT_CONNECTED) {
						aax := ax + rcGetDirOffsetX(2)
						aay := ay + rcGetDirOffsetY(2)
						aai := (chf cells[aax+aay*w] index as Int) + rcGetCon(asz, 2)
						if (src[aai]+3 < src[i])
							src[i] = src[aai]+3
					}
				}
			}
		}
	}
	
	// Pass 2
	for (y: Int = h-1; y >= 0; --y) {
		for (x: Int = w-1; x >= 0; --x) {
			c := chf cells[x+y*w]
			for (i: Int in (c index as Int)..(c index+c count) as Int) {
				s := chf spans[i]
				
				if (rcGetCon(s, 2) != RC_NOT_CONNECTED) {
					// (1,0)
					ax := x + rcGetDirOffsetX(2)
					ay := y + rcGetDirOffsetY(2)
					ai := (chf cells[ax+ay*w] index as Int) + rcGetCon(s, 2)
					asz := chf spans[ai]
					if (src[ai]+2 < src[i])
						src[i] = src[ai]+2
					
					// (1,1)
					if (rcGetCon(asz, 1) != RC_NOT_CONNECTED) {
						aax := ax + rcGetDirOffsetX(1)
						aay := ay + rcGetDirOffsetY(1)
						aai := (chf cells[aax+aay*w] index as Int) + rcGetCon(asz, 1)
						if (src[aai]+3 < src[i])
							src[i] = src[aai]+3
					}
				}
				if (rcGetCon(s, 1) != RC_NOT_CONNECTED) {
					// (0,1)
					ax := x + rcGetDirOffsetX(1)
					ay := y + rcGetDirOffsetY(1)
					ai := (chf cells[ax+ay*w] index as Int) + rcGetCon(s, 1)
					asz := chf spans[ai]
					if (src[ai]+2 < src[i])
						src[i] = src[ai]+2
					
					// (-1,1)
					if (rcGetCon(asz, 0) != RC_NOT_CONNECTED) {
						aax := ax + rcGetDirOffsetX(0)
						aay := ay + rcGetDirOffsetY(0)
						aai := (chf cells[aax+aay*w] index as Int) + rcGetCon(asz, 0)
						if (src[aai]+3 < src[i])
							src[i] = src[aai]+3
					}
				}
			}
		}
	}	
	
	maxDist = 0
	for (i: Int in 0..chf spanCount)
		maxDist = rcMax(src[i], maxDist)
	
	return src
}

boxBlur: static func(chf: RCCompactHeightfield@, thr: Int, src, dst: UShort*) -> UShort* {
	w := chf width
	h := chf height
	
	thr *= 2
	for (y: Int in 0..h) {
		for (x: Int in 0..w) {
			c := chf cells[x+y*w]
			for (i: Int in (c index as Int)..(c index + c count) as Int) {
				s := chf spans[i]
				cd := src[i] as Int
				if (cd <= thr) {
					dst[i] = cd
					continue
				}
				
				d := cd
				for (dir: Int in 0..4) {
					if (rcGetCon(s, dir) != RC_NOT_CONNECTED) {
						ax := x + rcGetDirOffsetX(dir)
						ay := y + rcGetDirOffsetY(dir)
						ai := (chf cells[ax+ay*w] index as Int) + rcGetCon(s, dir)
						d += src[ai] as Int
						
						asz := chf spans[ai]
						dir2 := (dir+1) & 0x3
						if (rcGetCon(asz, dir2) != RC_NOT_CONNECTED) {
							ax2 := ax + rcGetDirOffsetX(dir2)
							ay2 := ay + rcGetDirOffsetY(dir2)
							ai2 := (chf cells[ax2+ay2*w] index as Int) + rcGetCon(asz, dir2)
							d += (src[ai2] as Int)
						} else {
							d += cd
						}
					} else {
						d += cd*2
					}
				}
				dst[i] = ((d+5)/9) as UShort
			}
		}
	}
	return dst
}


floodRegion: static func(x, y, i: Int, level, r: UShort, chf: RCCompactHeightfield@, srcReg, srcDist: UShort*, stack: RCIntArray@) -> Bool{
	w := chf width
	area := chf areas[i]
	
	// Flood fill mark region.
	stack resize(0)
	stack push(x as Int)
	stack push(y as Int)
	stack push(i as Int)
	srcReg[i] = r
	srcDist[i] = 0
	
	lev: UShort = level >= 2 ? level-2 : 0
	count := 0
	
	while (stack size() > 0) {
		ci := stack pop()
		cy := stack pop()
		cx := stack pop()
		
		cs := chf spans[ci]
		
		// Check if any of the neighbours already have a valid region set.
		ar: UShort = 0
		for (dir: Int in 0..4) {
			// 8 connected
			if (rcGetCon(cs, dir) != RC_NOT_CONNECTED) {
				ax := cx + rcGetDirOffsetX(dir)
				ay := cy + rcGetDirOffsetY(dir)
				ai := (chf cells[ax+ay*w] index as Int) + rcGetCon(cs, dir)
				if (chf areas[ai] != area)
					continue
				nr := srcReg[ai]
				if (nr != 0 && nr != r)
					ar = nr
				
				asz := chf spans[ai]
				
				dir2 := (dir+1) & 0x3
				if (rcGetCon(asz, dir2) != RC_NOT_CONNECTED) {
					ax2 := ax + rcGetDirOffsetX(dir2)
					ay2 := ay + rcGetDirOffsetY(dir2)
					ai2 := (chf cells[ax2+ay2*w] index asz Int) + rcGetCon(asz, dir2)
					if (chf areas[ai2] != area)
						continue
					nr := srcReg[ai2]
					if (nr != 0 && nr != r)
						ar = nr
				}				
			}
		}
		if (ar != 0) {
			srcReg[ci] = 0
			continue
		}
		count+=1
		
		// Expand neighbours.
		for (dir: Int in 0..4) {
			if (rcGetCon(cs, dir) != RC_NOT_CONNECTED) {
				ax := cx + rcGetDirOffsetX(dir)
				ay := cy + rcGetDirOffsetY(dir)
				ai := (chf cells[ax+ay*w] index as Int) + rcGetCon(cs, dir)
				if (chf areas[ai] != area)
					continue
				if (chf dist[ai] >= lev) {
					if (srcReg[ai] == 0) {
						srcReg[ai] = r
						srcDist[ai] = 0
						stack push(ax)
						stack push(ay)
						stack push(ai)
					}
				}
			}
		}
	}
	
	return count > 0
}

expandRegions: static func(maxIter: Int, level: UShort, chf: RCCompactHeightfield@, srcReg, srcDist, dstReg, dstDist: UShort*, stack: RCIntArray@) -> UShort* {
	w := chf width
	h := chf height
	
	// Find cells revealed by the raised level.
	stack resize(0)
	for (y: Int in 0..h) {
		for (x: Int in 0..w) {
			c := chf cells[x+y*w]
			for (i: Int in (c index as Int)..(c index + c count) as Int) {
				if (chf dist[i] >= level && srcReg[i] == 0 && chf areas[i] != RC_NULL_AREA) {
					stack push(x)
					stack push(y)
					stack push(i)
				}
			}
		}
	}
	
	iter := 0
	while (stack size() > 0) {
		failed := 0
		memcpy(dstReg, srcReg, sizeof(UShort)*chf spanCount)
		memcpy(dstDist, srcDist, sizeof(UShort)*chf spanCount)
		
		for (j: Int = 0; j < stack size(); j += 3) {
			x := stack[j+0]
			y := stack[j+1]
			i := stack[j+2]
			if (i < 0) {
				failed+=1
				continue
			}
			
			r := srcReg[i]
			d2: UShort = 0xffff
			area := chf areas[i]
			s := chf spans[i]
			for (dir: Int in 0..4) {
				if (rcGetCon(s, dir) == RC_NOT_CONNECTED) continue
				ax := x + rcGetDirOffsetX(dir)
				ay := y + rcGetDirOffsetY(dir)
				ai := (chf cells[ax+ay*w] index as Int) + rcGetCon(s, dir)
				if (chf areas[ai] != area) continue
				if (srcReg[ai] > 0 && (srcReg[ai] & RC_BORDER_REG) == 0) {
					if ((srcDist[ai]+2 as Int) < (d2 as Int)) {
						r = srcReg[ai]
						d2 = srcDist[ai]+2
					}
				}
			}
			if (r) {
				stack[j+2] = -1// mark as used
				dstReg[i] = r
				dstDist[i] = d2
			} else {
				failed+=1
			}
		}
		
		// rcSwap source and dest.
		rcSwap(srcReg, dstReg)
		rcSwap(srcDist, dstDist)
		
		if (failed*3 == stack size())
			break
		
		if (level > 0) {
			iter+=1
			if (iter >= maxIter)
				break
		}
	}
	
	return srcReg
}

RCRegion: class {
	init: func() {}
	
	count: Int
	id: UShort
	area: UInt8
	remap: Bool
	connections: RCIntArray
	floors: RCIntArray
}

removeAdjacentNeighbours: static func(reg: RCRegion@) {
	// Remove adjacent duplicates.
	i := 0
	while (i < reg connections size() && reg connections size() > 1) {
		ni := (i+1) % reg connections size()
		if (reg connections[i] == reg connections[ni]) {
			// Remove duplicate
			for (j: Int in i..(reg connections size()-1))
				reg connections[j] = reg connections[j+1]
			reg connections pop()
		} else {
			i+=1
		}
	}
}

replaceNeighbour: static func(reg: RCRegion@, oldId, newId: UShort) {
	neiChanged := false
	for (i: Int in 0..(reg connections size())) {
		if (reg connections[i] == oldId) {
			reg connections[i] = newId
			neiChanged = true
		}
	}
	for (i: Int in 0..(reg floors size())) {
		if (reg floors[i] == oldId)
			reg floors[i] = newId
	}
	if (neiChanged)
		removeAdjacentNeighbours(reg)
}

canMergeWithRegion: static func(rega, regb: RCRegion@) -> Bool {
	if (rega area != regb area)
		return false
	n := 0
	for (i: Int in 0..(rega connections size())) {
		if (rega connections[i] == regb id)
			n+=1
	}
	if (n > 1)
		return false
	for (i: Int in 0..(rega floors size())) {
		if (rega floors[i] == regb id)
			return false
	}
	return true
}

addUniqueFloorRegion: static func(reg: RCRegion@, n: UShort) {
	for (i: Int in 0..(reg floors size()))
		if (reg floors[i] == n)
			return
	reg floors push(n)
}

mergeRegions: static func(rega, regb: RCRegion@) -> Bool {
	aid := rega id
	bid := regb id
	
	// Duplicate current neighbourhood.
	acon := RCIntArray new(rega connections size())
	for (i: Int in 0..(rega connections size()))
		acon[i] = rega connections[i]
	bcon := regb connections
	
	// Find insertion point on A.
	insa := -1
	for (i: Int in 0..(acon size())) {
		if (acon[i] == bid) {
			insa = i
			break
		}
	}
	if (insa == -1)
		return false
	
	// Find insertion point on B.
	insb := -1
	for (i: Int in 0..(bcon size())) {
		if (bcon[i] == aid) {
			insb = i
			break
		}
	}
	if (insb == -1)
		return false
	
	// Merge neighbours.
	rega connections resize(0)
	ni := acon size()
	for (i: Int in 0..(ni-1))
		rega connections push(acon[(insa+1+i) % ni])
		
	ni = bcon size()
	for (i: Int in 0..(ni-1))
		rega connections push(bcon[(insb+1+i) % ni])
	
	removeAdjacentNeighbours(rega)
	
	for (j: Int in 0..(regb floors size()))
		addUniqueFloorRegion(rega, regb floors[j])
	rega count += regb count
	regb count = 0
	regb connections resize(0)
	
	return true
}

isRegionConnectedToBorder: static func(reg: RCRegion@) -> Bool {
	// Region is connected to border if
	// one of the neighbours is null id.
	for (i: Int in 0..(reg connections size())) {
		if (reg connections[i] == 0)
			return true
	}
	return false
}

isSolidEdge: static func(chf: RCCompactHeightfield@, srcReg: UShort*, x, y, i, dir: Int) -> Bool {
	s := chf spans[i]
	r: UShort = 0
	if (rcGetCon(s, dir) != RC_NOT_CONNECTED) {
		ax := x + rcGetDirOffsetX(dir)
		ay := y + rcGetDirOffsetY(dir)
		ai := (chf cells[ax+ay*chf width] index as Int) + rcGetCon(s, dir)
		r = srcReg[ai]
	}
	if (r == srcReg[i])
		return false
	return true
}

walkContour: static func(x, y, i, dir: Int, chf: RCCompactHeightfield@, srcReg: UShort*, cont: RCIntArray@) {
	startDir := dir
	starti := i

	ss := chf spans[i]
	curReg: UShort = 0
	if (rcGetCon(ss, dir) != RC_NOT_CONNECTED) {
		ax := x + rcGetDirOffsetX(dir)
		ay := y + rcGetDirOffsetY(dir)
		ai := (chf cells[ax+ay*chf width] index as Int) + rcGetCon(ss, dir)
		curReg = srcReg[ai]
	}
	cont push(curReg)
	
	iter := 1 // 0
	while (iter < 40000) {
		s := chf spans[i]
		
		if (isSolidEdge(chf, srcReg, x, y, i, dir)) {
			// Choose the edge corner
			UShort r = 0
			if (rcGetCon(s, dir) != RC_NOT_CONNECTED) {
				ax := x + rcGetDirOffsetX(dir)
				ay := y + rcGetDirOffsetY(dir)
				ai := (chf cells[ax+ay*chf width] index as Int) + rcGetCon(s, dir)
				r = srcReg[ai]
			}
			if (r != curReg) {
				curReg = r
				cont push(curReg)
			}
			
			dir = (dir+1) & 0x3		// Rotate CW
		} else {
			ni := -1
			nx := x + rcGetDirOffsetX(dir)
			ny := y + rcGetDirOffsetY(dir)
			if (rcGetCon(s, dir) != RC_NOT_CONNECTED) {
				nc := chf cells[nx+ny*chf width]
				ni = (nc index as Int) + rcGetCon(s, dir)
			}
			if (ni == -1) {
				// Should not happen.
				return
			}
			x = nx
			y = ny
			i = ni
			dir = (dir+3) & 0x3		// Rotate CCW
		}
		
		if (starti == i && startDir == dir) {
			break
		}
		iter+=1
	}

	// Remove adjacent duplicates.
	if (cont size() > 1) {
		i := 0
		while (i < cont size()) {
			ni = (i+1) % (cont size())
			if (cont[i] == cont[ni]) {
				for (j: Int in i..(cont size()-1))
					cont[j] = cont[j+1]
				cont pop()
			} else
				i+=1
		}
	}
}

filterSmallRegions: static func(minRegionSize, mergeRegionSize: Int, maxRegionId: UShort@, chf: RCCompactHeightfield@, srcReg: UShort*) -> Bool {
	w := chf width
	h := chf height
	
	nreg := maxRegionId+1
	regions := rcAllocArray(RCRegion, nreg)
	if (!regions) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "filterSmallRegions: Out of memory 'regions' (%d).", nreg)
		return false
	}
	
	for (i: Int in 0..nreg)
		regions[i] id = i as UShort
	
	// Find edge of a region and find connections around the contour.
	for (y: Int in 0..h) {
		for (x: Int in 0..w) {
			c := chf cells[x+y*w]
			ni := (c index + c count) as Int
			for (i: Int in (c index as Int)..ni) {
				r := srcReg[i]
				if (r == 0 || r >= nreg)
					continue
				
				reg := regions[r]
				reg count+=1
				
				// Update floors.
				for (j: Int in (c index as Int)..ni) {
					if (i == j) continue
					floorId := srcReg[j]
					if (floorId == 0 || floorId >= nreg)
						continue
					addUniqueFloorRegion(reg, floorId)
				}
				
				// Have found contour
				if (reg connections size() > 0)
					continue
				
				reg area = chf areas[i]
				
				// Check if this cell is next to a border.
				ndir := -1
				for (dir: Int in 0..4) {
					if (isSolidEdge(chf, srcReg, x, y, i, dir)) {
						ndir = dir
						break
					}
				}
				
				if (ndir != -1) {
					// The cell is at border.
					// Walk around the contour to find all the neighbours.
					walkContour(x, y, i, ndir, chf, srcReg, reg connections)
				}
			}
		}
	}
	
	// Remove too small unconnected regions.
	for (i: Int in 0..nreg) {
		reg := regions[i]
		if (reg id == 0 || (reg id & RC_BORDER_REG))
			continue
		if (reg count == 0)
			continue
		
		if (reg connections size() == 1 && reg connections[0] == 0) {
			if (reg count < minRegionSize) {
				// Non-connected small region, remove.
				reg count = 0
				reg id = 0
			}
		}
	}
	
	// Merge too small regions to neighbour regions.
	mergeCount := 1 //0
	while (mergeCount > 0) { //do
		mergeCount = 0
		for (i: Int in 0..nreg) {
			reg := regions[i]
			if (reg id == 0 || (reg id & RC_BORDER_REG))
				continue
			if (reg count == 0)
				continue
			
			// Check to see if the region should be merged.
			if (reg count > mergeRegionSize && isRegionConnectedToBorder(reg))
				continue
			
			// Small region with more than 1 connection.
			// Or region which is not connected to a border at all.
			// Find smallest neighbour region that connects to this one.
			smallest := 0xfffffff
			mergeId := reg id
			for (j: Int in 0..(reg connections size())) {
				if (reg connections[j] & RC_BORDER_REG) continue
				mreg := regions[reg connections[j]]
				if (mreg id == 0 || (mreg id & RC_BORDER_REG)) continue
				if (mreg count < smallest &&
					canMergeWithRegion(reg, mreg) &&
					canMergeWithRegion(mreg, reg)) {
					smallest = mreg count
					mergeId = mreg id
				}
			}
			// Found new id.
			if (mergeId != reg id) {
				oldId := reg id
				target := regions[mergeId]
				
				// Merge neighbours.
				if (mergeRegions(target, reg)) {
					// Fixup regions pointing to current region.
					for (j: Int in 0..nreg) {
						if (regions[j] id == 0 || (regions[j] id & RC_BORDER_REG)) continue
						// If another region was already merged into current region
						// change the nid of the previous region too.
						if (regions[j] id == oldId)
							regions[j] id = mergeId
						// Replace the current region with the new one if the
						// current regions is neighbour.
						replaceNeighbour(regions[j], oldId, mergeId)
					}
					mergeCount+=1
				}
			}
		}
	} // while (mergeCount > 0)
	
	// Compress region Ids.
	for (i: Int in 0..nreg) {
		regions[i] remap = false
		if (regions[i] id == 0) continue// Skip nil regions.
		if (regions[i] id & RC_BORDER_REG) continue// Skip external regions.
		regions[i] remap = true
	}
	
	regIdGen: UShort = 0
	for (i: Int in 0..nreg) {
		if (!regions[i] remap)
			continue
		oldId := regions[i] id
		regIdGen+=1
		newId := regIdGen
		for (j: Int in i..nreg) {
			if (regions[j] id == oldId) {
				regions[j] id = newId
				regions[j] remap = false
			}
		}
	}
	maxRegionId = regIdGen
	
	// Remap regions.
	for (i: Int in 0..(chf spanCount)) {
		if ((srcReg[i] & RC_BORDER_REG) == 0)
			srcReg[i] = regions[srcReg[i]] id
	}
	delete [] regions
	return true
}

rcBuildDistanceField: func(chf: RCCompactHeightfield@) -> Bool{
	startTime := rcGetPerformanceTimer()
	
	if (chf dist) {
		delete [] chf dist
		chf dist = 0
	}
	
	dist0 := rcAllocArray(UShort, chf spanCount)
	if (!dist0) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildDistanceField: Out of memory 'dist0' (%d).", chf spanCount)
		return false
	}
	dist1 := rcAllocArray(UShort, chf spanCount)
	if (!dist1) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildDistanceField: Out of memory 'dist1' (%d).", chf spanCount)
		delete [] dist0
		return false
	}
	
	src := dist0
	dst := dist1
	maxDist: UShort = 0
	
	diststartTime := rcGetPerformanceTimer()
	if (calculateDistanceField(chf, src, dst, maxDist) != src)
		rcSwap(src, dst)
	
	chf maxDistance = maxDist
	distendTime := rcGetPerformanceTimer()
	
	blurstartTime := rcGetPerformanceTimer()
	// Blur
	if (boxBlur(chf, 1, src, dst) != src)
		rcSwap(src, dst)
	
	// Store distance.
	chf dist = src
	blurendTime := rcGetPerformanceTimer()
	delete [] dst
	
	endTime := rcGetPerformanceTimer()
	
/*	if (rcGetLog()) {
		rcGetLog() log(RCLogCategory RC_LOG_PROGRESS, "Build distance field: %.3f ms", rcGetDeltaTimeUsec(startTime, endTime)/1000.0f)
		rcGetLog() log(RCLogCategory RC_LOG_PROGRESS, " - dist: %.3f ms", rcGetDeltaTimeUsec(distStartTime, distEndTime)/1000.0f)
		rcGetLog() log(RCLogCategory RC_LOG_PROGRESS, " - blur: %.3f ms", rcGetDeltaTimeUsec(blurStartTime, blurEndTime)/1000.0f)
	}*/
	if (rcGetBuildTimes()) {
		rcGetBuildTimes() buildDistanceField += rcGetDeltaTimeUsec(startTime, endTime)
		rcGetBuildTimes() buildDistanceFieldDist += rcGetDeltaTimeUsec(distStartTime, distEndTime)
		rcGetBuildTimes() buildDistanceFieldBlur += rcGetDeltaTimeUsec(blurStartTime, blurEndTime)
	}
	
	return true
}

paintRectRegion: static func(minx, maxx, miny, maxy: Int, regId: UShort, chf: RCCompactHeightfield@, srcReg: UShort*) {
	w = chf width
	for (y: Int in miny..maxy) {
		for (x: Int in minx..maxx) {
			c := chf cells[x+y*w]
			for (i: Int in (c index as Int)..(c index+c count) as Int) {
				if (chf areas[i] != RC_NULL_AREA)
					srcReg[i] = regId
			}
		}
	}
}


RC_NULL_NEI: static UShort = 0xffff

RCSweepSpan: class {
	rid: UShort		// row id
	id: UShort		// region id
	ns: UShort		// number samples
	nei: UShort		// neighbour id
}

rcBuildRegionsMonotone: func(chf: RCCompactHeightfield@, borderSize, minRegionSize, mergeRegionSize: Int) -> Bool {
	startTime := rcGetPerformanceTimer()
	
	w := chf width
	h := chf height
	id: UShort = 1
	
	if (chf regs) {
		delete [] chf regs
		chf regs = 0
	}
	
	srcReg := rcAllocArray(UShort, chf spanCount)
	if (!srcReg) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildRegionsMonotone: Out of memory 'src' (%d).", chf spanCount)
		return false
	}
	memset(srcReg, 0, sizeof(UShort)*chf spanCount)
	
	sweeps := rcAllocArray(RCSweepSpan, rcMax(chf width, chf height))
	if (!sweeps) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildRegionsMonotone: Out of memory 'sweeps' (%d).", chf width)
		return false
	}
	
	// Mark border regions.
	if (borderSize) {
		paintRectRegion(0, borderSize, 0, h, id|RC_BORDER_REG, chf, srcReg); id+=1
		paintRectRegion(w-borderSize, w, 0, h, id|RC_BORDER_REG, chf, srcReg); id+=1
		paintRectRegion(0, w, 0, borderSize, id|RC_BORDER_REG, chf, srcReg); id+=1
		paintRectRegion(0, w, h-borderSize, h, id|RC_BORDER_REG, chf, srcReg); id+=1
	}
	
	prev := RCIntArray new(256)
	
	// Sweep one line at a time.
	for (y: Int in borderSize..(h-borderSize)) {
		// Collect spans from this row.
		prev resize(id+1)
		memset(prev[0]&, 0, sizeof(Int)*id)
		rid: UShort = 1
		
		for (x: Int in borderSize..(w-borderSize)) {
			c := chf cells[x+y*w]
			for (i: Int in (c index as Int)..(c index+c count) as Int) {
				s := chf spans[i]
				if (chf areas[i] == RC_NULL_AREA) continue
				
				// -x
				previd: UShort = 0
				if (rcGetCon(s, 0) != RC_NOT_CONNECTED) {
					ax := x + rcGetDirOffsetX(0)
					ay := y + rcGetDirOffsetY(0)
					ai := (chf cells[ax+ay*w] index as Int) + rcGetCon(s, 0)
					if ((srcReg[ai] & RC_BORDER_REG) == 0 && chf areas[i] == chf areas[ai])
						previd = srcReg[ai]
				}
				
				if (!previd) {
					previd = rid
					rid+=1
					sweeps[previd] rid = previd
					sweeps[previd] ns = 0
					sweeps[previd] nei = 0
				}

				// -y
				if (rcGetCon(s,3) != RC_NOT_CONNECTED) {
					ax := x + rcGetDirOffsetX(3)
					ay := y + rcGetDirOffsetY(3)
					ai := (chf cells[ax+ay*w] index as Int) + rcGetCon(s, 3)
					if (srcReg[ai] && (srcReg[ai] & RC_BORDER_REG) == 0 && chf areas[i] == chf areas[ai]) {
						nr := srcReg[ai]
						if (!sweeps[previd] nei || sweeps[previd] nei == nr) {
							sweeps[previd] nei = nr
							sweeps[previd] ns+=1
							prev[nr]+=1
						} else {
							sweeps[previd] nei = RC_NULL_NEI
						}
					}
				}
				srcReg[i] = previd
			}
		}
		
		// Create unique ID.
		for (i: Int in 1..rid) {
			if (sweeps[i] nei != RC_NULL_NEI && sweeps[i] nei != 0 && prev[sweeps[i] nei] == (sweeps[i] ns as Int)) {
				sweeps[i] id = sweeps[i] nei
			} else {
				sweeps[i] id = id
				id+=1
			}
		}
		
		// Remap IDs
		for (x: Int in borderSize..(w-borderSize)) {
			c := chf cells[x+y*w]
			for (i: Int in (c index as Int)..(c index+c count) as Int) {
				if (srcReg[i] > 0 && srcReg[i] < rid)
					srcReg[i] = sweeps[srcReg[i]] id
			}
		}
	}
	
	filterstartTime := rcGetPerformanceTimer()
	// Filter out small regions.
	chf maxRegions = id
	if (!filterSmallRegions(minRegionSize, mergeRegionSize, chf maxRegions, chf, srcReg))
		return false
	
	filterendTime := rcGetPerformanceTimer()
	
	// Store the result out.
	chf regs = srcReg
	srcReg = 0
	
	endTime := rcGetPerformanceTimer()
	if (rcGetBuildTimes()) {
		rcGetBuildTimes() buildRegions += rcGetDeltaTimeUsec(startTime, endTime)
		rcGetBuildTimes() buildRegionsFilter += rcGetDeltaTimeUsec(filterStartTime, filterEndTime)
	}

	return true
}

rcBuildRegions: func(chf: RCCompactHeightfield@, borderSize, minRegionSize, mergeRegionSize: Int) -> Bool {
	startTime := rcGetPerformanceTimer()
	w := chf width
	h := chf height
	
	if (!chf regs) {
		chf regs = rcAllocArray(UShort, chf spanCount)
		if (!chf regs) {
			if (rcGetLog())
				rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildRegions: Out of memory 'chf.reg' (%d).", chf spanCount)
			return false
		}
	}
	
	tmp := rcAllocArray(UShort, chf spanCount*4)
	if (!tmp) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildRegions: Out of memory 'tmp' (%d).", chf spanCount*4)
		return false
	}
	
	regstartTime := rcGetPerformanceTimer()
	
	stack := RCIntArray new(1024)
	visited := RCIntArray new(1024)
	
	srcReg: UShort* = tmp
	srcDist: UShort* = tmp+chf spanCount
	dstReg: UShort* = tmp+chf spanCount*2
	dstDist: UShort* = tmp+chf spanCount*3
	
	memset(srcReg, 0, sizeof(UShort)*chf spanCount)
	memset(srcDist, 0, sizeof(UShort)*chf spanCount)
	
	regionId: UShort = 1
	level: UShort = (chf maxDistance+1) & ~1
	
	// TODO: Figure better formula, expandIters defines how much the 
	// watershed "overflows" and simplifies the regions. Tying it to
	// agent radius was usually good indication how greedy it could be.
//	Int expandIters = 4 + walkableRadius * 2
	expandIters := 8
	
	// Mark border regions.
	paintRectRegion(0, borderSize, 0, h, regionId|RC_BORDER_REG, chf, srcReg); regionId+=1
	paintRectRegion(w-borderSize, w, 0, h, regionId|RC_BORDER_REG, chf, srcReg); regionId+=1
	paintRectRegion(0, w, 0, borderSize, regionId|RC_BORDER_REG, chf, srcReg); regionId+=1
	paintRectRegion(0, w, h-borderSize, h, regionId|RC_BORDER_REG, chf, srcReg); regionId+=1
	
	expTime = 0
	floodTime = 0
	
	while (level > 0) {
		level = level >= 2 ? level-2 : 0
		expstartTime := rcGetPerformanceTimer()
		
		// Expand current regions until no empty connected cells found.
		if (expandRegions(expandIters, level, chf, srcReg, srcDist, dstReg, dstDist, stack) != srcReg) {
			rcSwap(srcReg, dstReg)
			rcSwap(srcDist, dstDist)
		}
		expTime += rcGetPerformanceTimer() - expStartTime
		
		floodstartTime := rcGetPerformanceTimer()
		// Mark new regions with IDs.
		for (y: Int in 0..h) {
			for (x: Int in 0..w) {
				c := chf cells[x+y*w]
				for (i: Int in (c index as Int)..(c index+c count) as Int) {
					if (chf dist[i] < level || srcReg[i] != 0 || chf areas[i] == RC_NULL_AREA)
						continue
					
					if (floodRegion(x, y, i, level, regionId, chf, srcReg, srcDist, stack))
						regionId+=1
				}
			}
		}
		floodTime += rcGetPerformanceTimer() - floodStartTime
	}
	
	// Expand current regions until no empty connected cells found.
	if (expandRegions(expandIters*8, 0, chf, srcReg, srcDist, dstReg, dstDist, stack) != srcReg) {
		rcSwap(srcReg, dstReg)
		rcSwap(srcDist, dstDist)
	}
	regendTime := rcGetPerformanceTimer()
	
	filterstartTime := rcGetPerformanceTimer()
	// Filter out small regions.
	chf maxRegions = regionId
	if (!filterSmallRegions(minRegionSize, mergeRegionSize, chf maxRegions, chf, srcReg))
		return false
	
	filterendTime := rcGetPerformanceTimer()
		
	// Write the result out.
	memcpy(chf regs, srcReg, sizeof(UShort)*chf spanCount)
	endTime := rcGetPerformanceTimer()
	
/*	if (rcGetLog()) {
		rcGetLog() log(RCLogCategory RC_LOG_PROGRESS, "Build regions: %.3f ms", rcGetDeltaTimeUsec(startTime, endTime)/1000.0f)
		rcGetLog() log(RCLogCategory RC_LOG_PROGRESS, " - reg: %.3f ms", rcGetDeltaTimeUsec(regStartTime, regEndTime)/1000.0f)
		rcGetLog() log(RCLogCategory RC_LOG_PROGRESS, " - exp: %.3f ms", rcGetDeltaTimeUsec(0, expTime)/1000.0f)
		rcGetLog() log(RCLogCategory RC_LOG_PROGRESS, " - flood: %.3f ms", rcGetDeltaTimeUsec(0, floodTime)/1000.0f)
		rcGetLog() log(RCLogCategory RC_LOG_PROGRESS, " - filter: %.3f ms", rcGetDeltaTimeUsec(filterStartTime, filterEndTime)/1000.0f)
	}
*/
	if (rcGetBuildTimes()) {
		rcGetBuildTimes() buildRegions += rcGetDeltaTimeUsec(startTime, endTime)
		rcGetBuildTimes() buildRegionsReg += rcGetDeltaTimeUsec(regStartTime, regEndTime)
		rcGetBuildTimes() buildRegionsExp += rcGetDeltaTimeUsec(0, expTime)
		rcGetBuildTimes() buildRegionsFlood += rcGetDeltaTimeUsec(0, floodTime)
		rcGetBuildTimes() buildRegionsFilter += rcGetDeltaTimeUsec(filterStartTime, filterEndTime)
	}
		
	return true
}

