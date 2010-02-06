use retour
import Retour/[Retour, Log, Timer]

getCornerHeight: static func(x, y, i, dir: Int, chf: RCCompactHeightfield@, isBorderVertex: Bool@) -> Int {
	s := chf spans[i]
	ch := s y as Int
	dirp := (dir + 1) & 0x3
	
	regs: UInt[4] = [0, 0, 0, 0]
	
	// Combine region and area codes in order to prevent
	// border vertices which are in between two areas to be removed. 
	regs[0] = chf regs[i] | (chf areas[i] << 16)
	
	if (rcGetCon(s, dir) != RC_NOT_CONNECTED) {
		ax := x + rcGetDirOffsetX(dir)
		ay := y + rcGetDirOffsetY(dir)
		ai = (chf cells[ax+ay*chf width] index as Int) + rcGetCon(s, dir)
		asp := chf spans[ai]
		ch = rcMax(ch, asp y as Int)
		regs[1] = chf regs[ai] | (chf areas[ai] << 16)
		if (rcGetCon(asp, dirp) != RC_NOT_CONNECTED) {
			ax2 := ax + rcGetDirOffsetX(dirp)
			ay2 := ay + rcGetDirOffsetY(dirp)
			ai2 := (chf cells[ax2+ay2*chf width] index as Int) + rcGetCon(asp, dirp)
			as2 := chf spans[ai2]
			ch = rcMax(ch, as2 y as Int)
			regs[2] = chf regs[ai2] | (chf areas[ai2] << 16)
		}
	}
	if (rcGetCon(s, dirp) != RC_NOT_CONNECTED) {
		ax = x + rcGetDirOffsetX(dirp)
		ay = y + rcGetDirOffsetY(dirp)
		ai = (chf cells[ax+ay*chf width] index as Int) + rcGetCon(s, dirp)
		asp := chf spans[ai]
		ch = rcMax(ch, asp y as Int)
		regs[3] = chf regs[ai] | (chf areas[ai] << 16)
		if (rcGetCon(asp, dir) != RC_NOT_CONNECTED) {
			ax2 := ax + rcGetDirOffsetX(dir)
			ay2 := ay + rcGetDirOffsetY(dir)
			ai2 := (chf cells[ax2+ay2*chf width] index as Int) + rcGetCon(asp, dir)
			as2 := chf spans[ai2]
			ch = rcMax(ch, as2 y as Int)
			regs[2] = chf regs[ai2] | (chf areas[ai2] << 16)
		}
	}
	
	// Check if the vertex is special edge vertex, these vertices will be removed later.
	for (j in 0..4) {
		a := j
		b := (j+1) & 0x3
		c := (j+2) & 0x3
		d := (j+3) & 0x3
		
		// The vertex is a border vertex there are two same exterior cells in a row,
		// followed by two interior cells and none of the regions are out of bounds.
		twoSameExts: Bool = (regs[a] & regs[b] & RC_BORDER_REG) != 0 && regs[a] == regs[b]
		twoInts: Bool = ((regs[c] | regs[d]) & RC_BORDER_REG) == 0
		intsSameArea: Bool = (regs[c]>>16) == (regs[d]>>16)
		noZeros: Bool = regs[a] != 0 && regs[b] != 0 && regs[c] != 0 && regs[d] != 0
		if (twoSameExts && twoInts && intsSameArea && noZeros) {
			isBorderVertex = true
			break
		}
	}
	
	return ch
}

walkContour: static func(x, y, i: Int, chf: RCCompactHeightfield@, flags: UInt8*, points: RCIntArray@) {
	// Choose the first non-connected edge
	dir: UInt8 = 0
	while ((flags[i] & (1 << dir)) == 0)
		dir += 1
	
	startDir := dir
	starti := i
	
	area := chf areas[i]
	
	iter := 1 // iter := 0
	while (iter < 40000) { //while (++iter < 40000)
		if (flags[i] & (1 << dir)) {
			// Choose the edge corner
			isBorderVertex := false
			isAreaBorder := false
			px := x
			py := getCornerHeight(x, y, i, dir, chf, isBorderVertex)
			pz := y
			match dir {
				case 0 => pz+=1; break
				case 1 => px+=1; pz+=1; break
				case 2 => px+=1; break
			}
			r := 0
			s := chf spans[i]
			if (rcGetCon(s, dir) != RC_NOT_CONNECTED) {
				ax := x + rcGetDirOffsetX(dir)
				ay := y + rcGetDirOffsetY(dir)
				ai := (chf cells[ax+ay*chf width] as Int) index + rcGetCon(s, dir)
				r = (int)chf regs[ai]
				if (area != chf areas[ai])
					isAreaBorder = true
			}
			if (isBorderVertex)
				r |= RC_BORDER_VERTEX
			if (isAreaBorder)
				r |= RC_AREA_BORDER
			points push(px)
			points push(py)
			points push(pz)
			points push(r)
			
			flags[i] &= ~(1 << dir); // Remove visited edges
			dir = (dir+1) & 0x3;  // Rotate CW
		} else {
			ni := -1
			nx := x + rcGetDirOffsetX(dir)
			ny := y + rcGetDirOffsetY(dir)
			s := chf spans[i]
			if (rcGetCon(s, dir) != RC_NOT_CONNECTED) {
				nc := chf cells[nx+ny*chf width]
				ni = nc index as Int + rcGetCon(s, dir)
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
		iter += 1
	}
}

distancePtSeg: static func(x, y, z, px, py, pz, qx, qy, qz: Int) -> Float {
/*	float pqx = (float)(qx - px)
	float pqy = (float)(qy - py)
	float pqz = (float)(qz - pz)
	float dx = (float)(x - px)
	float dy = (float)(y - py)
	float dz = (float)(z - pz)
	float d = pqx*pqx + pqy*pqy + pqz*pqz
	float t = pqx*dx + pqy*dy + pqz*dz
	if (d > 0)
		t /= d
	if (t < 0)
		t = 0
	else if (t > 1)
		t = 1
	
	dx = px + t*pqx - x
	dy = py + t*pqy - y
	dz = pz + t*pqz - z
	
	return dx*dx + dy*dy + dz*dz;*/

	pqx = (qx - px) as Float
	pqz = (qz - pz) as Float
	dx = (x - px) as Float
	dz = (z - pz) as Float
	d = pqx*pqx + pqz*pqz
	t = pqx*dx + pqz*dz
	if (d > 0)
		t /= d
	if (t < 0)
		t = 0
	else if (t > 1)
		t = 1
	
	dx = px + t*pqx - x
	dz = pz + t*pqz - z
	
	return dx*dx + dz*dz
}

simplifyContour: static func(points, simplified: RCIntArray@, maxError: Float, maxEdgeLen: Int) {
	// Add initial points.
	noConnections := true
	for (i := 0; i < (points size()); i+=4) {
		if ((points[i+3] & RC_CONTOUR_REG_MASK) != 0) {
			noConnections = false
			break
		}
	}
	
	if (noConnections) {
		// If there is no connections at all,
		// create some initial points for the simplification process.
		// Find lower-left and upper-right vertices of the contour.
		llx := points[0]
		lly := points[1]
		llz := points[2]
		lli := 0
		urx := points[0]
		ury := points[1]
		urz := points[2]
		uri := 0
		for (i := 0; i < points size(); i+=4) {
			int x = points[i+0]
			int y = points[i+1]
			int z = points[i+2]
			if (x < llx || (x == llx && z < llz)) {
				llx = x
				lly = y
				llz = z
				lli = i / 4
			}
			if (x >= urx || (x == urx && z > urz)) {
				urx = x
				ury = y
				urz = z
				uri = i / 4
			}
		}
		simplified push(llx)
		simplified push(lly)
		simplified push(llz)
		simplified push(lli)
		
		simplified push(urx)
		simplified push(ury)
		simplified push(urz)
		simplified push(uri)
	} else {
		// The contour has some portals to other regions.
		// Add a new point to every location where the region changes.
		for (i in 0..(points size()/4)) {
			ii := (i+1) % ni
			differentRegs: Bool = (points[i*4+3] & RC_CONTOUR_REG_MASK) != (points[ii*4+3] & RC_CONTOUR_REG_MASK)
			areaBorders: Bool = (points[i*4+3] & RC_AREA_BORDER) != (points[ii*4+3] & RC_AREA_BORDER)
			if (differentRegs || areaBorders) {
				simplified push(points[i*4+0])
				simplified push(points[i*4+1])
				simplified push(points[i*4+2])
				simplified push(i)
			}
		}
	}
	
	// Add points until all raw points are within
	// error tolerance to the simplified shape.
	pn := points size()/4
	for (i in 0..(simplified size()/4)) {
		ii := (i+1) % (simplified size()/4)
		
		ax := simplified[i*4+0]
		ay := simplified[i*4+1]
		az := simplified[i*4+2]
		ai := simplified[i*4+3]
		
		bx := simplified[ii*4+0]
		by := simplified[ii*4+1]
		bz := simplified[ii*4+2]
		bi := simplified[ii*4+3]

		// Find maximum deviation from the segment.
		maxd := 0.0
		maxi := -1
		ci, cinc, endi: Int
		
		// Traverse the segment in lexilogical order so that the
		// max deviation is calculated similarly when traversing
		// opposite segments.
		if (bx > ax || (bx == ax && bz > az)) {
			cinc = 1
			ci = (ai+cinc) % pn
			endi = bi
		} else {
			cinc = pn-1
			ci = (bi+cinc) % pn
			endi = ai
		}
		
		// Tesselate only outer edges oredges between areas.
		if ((points[ci*4+3] & RC_CONTOUR_REG_MASK) == 0 ||
			(points[ci*4+3] & RC_AREA_BORDER)) {
			while (ci != endi) {
				d := distancePtSeg(points[ci*4+0], points[ci*4+1]/4, points[ci*4+2],
										ax, ay/4, az, bx, by/4, bz)
				if (d > maxd) {
					maxd = d
					maxi = ci
				}
				ci = (ci+cinc) % pn
			}
		}
		
		// If the max deviation is larger than accepted error,
		// add new point, else continue to next segment.
		if (maxi != -1 && maxd > (maxError*maxError)) {
			// Add space for the new point.
			simplified resize(simplified size()+4)
			int n = simplified size()/4
			for (j: Int = n-1; j > i; --j) {
				simplified[j*4+0] = simplified[(j-1)*4+0]
				simplified[j*4+1] = simplified[(j-1)*4+1]
				simplified[j*4+2] = simplified[(j-1)*4+2]
				simplified[j*4+3] = simplified[(j-1)*4+3]
			}
			// Add the point.
			simplified[(i+1)*4+0] = points[maxi*4+0]
			simplified[(i+1)*4+1] = points[maxi*4+1]
			simplified[(i+1)*4+2] = points[maxi*4+2]
			simplified[(i+1)*4+3] = maxi
		} else {
			i += 1
		}
	}
	
	// Split too long edges.
	if (maxEdgeLen > 0) {
		for (i in 0..(simplified size()/4) ) {
			ii := (i+1) % (simplified size()/4)
			
			ax := simplified[i*4+0]
			az := simplified[i*4+2]
			ai := simplified[i*4+3]
			
			bx := simplified[ii*4+0]
			bz := simplified[ii*4+2]
			bi := simplified[ii*4+3]
			
			// Find maximum deviation from the segment.
			maxi := -1
			ci := (ai+1) % pn
			
			// Tesselate only outer edges.
			if ((points[ci*4+3] & RC_CONTOUR_REG_MASK) == 0) {
				dx := bx - ax
				dz := bz - az
				if (dx*dx + dz*dz > maxEdgeLen*maxEdgeLen) {
					n := bi < ai ? (bi+pn - ai) : (bi - ai)
					maxi = (ai + n/2) % pn
				}
			}
			
			// If the max deviation is larger than accepted error,
			// add new point, else continue to next segment.
			if (maxi != -1) {
				// Add space for the new point.
				simplified resize(simplified size()+4)
				n := (simplified size()/4)
				for (j := n-1; j > i; --j) {
					simplified[j*4+0] = simplified[(j-1)*4+0]
					simplified[j*4+1] = simplified[(j-1)*4+1]
					simplified[j*4+2] = simplified[(j-1)*4+2]
					simplified[j*4+3] = simplified[(j-1)*4+3]
				}
				// Add the point.
				simplified[(i+1)*4+0] = points[maxi*4+0]
				simplified[(i+1)*4+1] = points[maxi*4+1]
				simplified[(i+1)*4+2] = points[maxi*4+2]
				simplified[(i+1)*4+3] = maxi
			} else {
				i += 1
			}
		}
	}
	
	for (i: Int in 0..(simplified size()/4)) {
		// The edge vertex flag is taken from the current raw point,
		// and the neighbour region is take from the next raw point.
		ai := (simplified[i*4+3]+1) % pn
		bi := simplified[i*4+3]
		simplified[i*4+3] = (points[ai*4+3] & RC_CONTOUR_REG_MASK) | (points[bi*4+3] & RC_BORDER_VERTEX)
	}
	
}

removeDegenerateSegments: static func(simplified: RCIntArray@) {
	// Remove adjacent vertices which are equal on xz-plane,
	// or else the triangulator will get confused.
	for (i: Int in 0..(simplified size()/4)) {
		ni := i+1
		if (ni >= (simplified size()/4))
			ni = 0
			
		if (simplified[i*4+0] == simplified[ni*4+0] &&
			simplified[i*4+2] == simplified[ni*4+2]) {
			// Degenerate segment, remove.
			for (j: Int in i..(simplified size()/4-1)) {
				simplified[j*4+0] = simplified[(j+1)*4+0]
				simplified[j*4+1] = simplified[(j+1)*4+1]
				simplified[j*4+2] = simplified[(j+1)*4+2]
				simplified[j*4+3] = simplified[(j+1)*4+3]
			}
			simplified resize(simplified size()-4)
		}
	}
}

calcAreaOfPolygon2D: static func(verts: Int*, nverts: Int) -> Int {
	area := 0
	j := nverts-1
	for (i: Int in 0..nverts) {
		vi: Int* = verts[i*4]&
		vj: Int* = verts[j*4]&
		area += vi[0] * vj[2] - vj[0] * vi[2]
		j = i
	}
	return (area+1) / 2
}

getClosestIndices: static func(vertsa: Int*, nvertsa: Int, vertsb: Int*, nvertsb: Int, ia, ib: Int@) {
	closestDist := 0xfffffff
	for (i in 0..nvertsa) {
		va: Int* = vertsa[i*4]&
		for (j: Int in 0..nvertsb) {
			vb: Int* = vertsb[j*4]&
			dx := vb[0] - va[0]
			dz := vb[2] - va[2]
			d := dx*dx + dz*dz
			if (d < closestDist) {
				ia = i
				ib = j
				closestDist = d
			}
		}
	}
}

mergeContours: static func(ca, cb: RCContour@, ia, ib: Int) -> Bool {
	maxVerts := ca nverts + cb nverts + 2
	verts := Array<Int> new(maxVerts*4)
	if (!verts)
		return false
	
	nv := 0
	// Copy contour A.
	for (i := 0; i <= (ca nverts); i += 1) {
		dst: Int* = verts[nv*4]&
		src: Int* = ca verts[((ia+i)%ca nverts)*4]&
		dst[0] = src[0]
		dst[1] = src[1]
		dst[2] = src[2]
		dst[3] = src[3]
		nv += 1
	}
	
	// Copy contour B
	for (i := 0; i <= cb nverts; i += 1) {
		dst: Int* = verts[nv*4]&
		src: Int* = cb verts[((ib+i)%cb nverts)*4]&
		dst[0] = src[0]
		dst[1] = src[1]
		dst[2] = src[2]
		dst[3] = src[3]
		nv += 1
	}
	
	//delete [] ca verts
	ca verts = verts
	ca nverts = nv
	
	//delete [] cb verts
	cb verts = null
	cb nverts = null
	
	return true
}

rcBuildContours: func(chf: RCCompactHeightfield@, maxError: Float, maxEdgeLen: Int, cset: RCContourSet@) -> Bool {
	w := chf width
	h := chf height
	
	startTime := rcGetPerformanceTimer()
	
	vcopy(cset bmin, chf bmin)
	vcopy(cset bmax, chf bmax)
	cset cs = chf cs
	cset ch = chf ch
	
	maxContours := chf maxRegions*2
	cset conts = new RCContour[maxContours]
	if (!cset conts)
		return false
	cset nconts = 0
	
	flags := rcAllocArray(UInt8, chf spanCount)
	if (!flags) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildContours: Out of memory 'flags'.")
		return false
	}
	traceStartTime := rcGetPerformanceTimer()
	
	// Mark boundaries.
	for (y in 0..h) {
		for (x in 0..w) {
			c := chf cells[x+y*w]
			for (i: Int in (c index as Int)..(c index + c count as Int)) {
				res: UInt8 = 0
				s := chf spans[i]
				if (!chf regs[i] || (chf regs[i] & RC_BORDER_REG)) {
					flags[i] = 0
					continue
				}
				for (dir: Int in 0..4) {
					r: UShort = 0
					if (rcGetCon(s, dir) != RC_NOT_CONNECTED) {
						ax := x + rcGetDirOffsetX(dir)
						ay := y + rcGetDirOffsetY(dir)
						ai := (chf cells[ax+ay*w] index as Int) + rcGetCon(s, dir)
						r = chf regs[ai]
					}
					if (r == chf regs[i])
						res |= (1 << dir)
				}
				flags[i] = res ^ 0xf	// Inverse, mark non connected edges.
			}
		}
	}
	traceEndTime := rcGetPerformanceTimer()
	
	simplifyStartTime := rcGetPerformanceTimer()
	verts := RCIntArray new(256)
	simplified := RCIntArray new(64)
	
	for (y in 0..h) {
		for (x in 0..w) {
			c := chf cells[x+y*w]
			for (i in (c index as Int)..(c index + c count as Int)) {
				if (flags[i] == 0 || flags[i] == 0xf) {
					flags[i] = 0
					continue
				}
				reg := chf regs[i]
				if (!reg || (reg & RC_BORDER_REG))
					continue
				area := chf areas[i]
				
				verts resize(0)
				simplified resize(0)
				walkContour(x, y, i, chf, flags, verts) 
				simplifyContour(verts, simplified, maxError, maxEdgeLen)
				removeDegenerateSegments(simplified)
				
				// Store region contour remap info.
				// Create contour.
				if (simplified size()/4 >= 3) {
					if (cset nconts > maxContours) {
						if (rcGetLog())
							rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildContours: Too many contours %d, max %d.", cset nconts, maxContours)
						return false
					}
					
					cont := cset conts[cset nconts]&
					cset nconts += 1
					
					cont nverts = simplified size()/4
					cont verts = Array<Int> new(cont nverts*4)
					memcpy(cont verts, simplified[0]&, sizeof(Int)*cont nverts*4)
					
					cont nrverts = verts size()/4
					cont rverts = Array<Int> new(cont nrverts*4)
					memcpy(cont rverts, verts[0]&, sizeof(Int)*cont nrverts*4)
					
/*					cont cx = cont cy = cont cz = 0
					for (i := 0; i < cont nverts; ++i)
					{
						cont cx += cont verts[i*4+0]
						cont cy += cont verts[i*4+1]
						cont cz += cont verts[i*4+2]
					}
					cont cx /= cont nverts
					cont cy /= cont nverts
					cont cz /= cont nverts;*/
					
					cont reg = reg
					cont area = area
				}
			}
		}
	}
	
	// Check and merge droppings.
	// Sometimes the previous algorithms can fail and create several countours
	// per area. This pass will try to merge the holes into the main region.
	for (i in 0..cset nconts) {
		cont := cset conts[i]
		// Check if the contour is would backwards.
		if (calcAreaOfPolygon2D(cont verts, cont nverts) < 0) {
			// Find another contour which has the same region ID.
			mergeIdx := -1
			for (j in 0..(cset nconts)) {
				if (i == j) continue
				if (cset conts[j] nverts && (cset conts[j] reg) == (cont reg)) {
					// Make sure the polygon is correctly oriented.
					if (calcAreaOfPolygon2D(cset conts[j] verts, cset conts[j] nverts)) {
						mergeIdx = j
						break
					}
				}
			}
			if (mergeIdx == -1) {
				if (rcGetLog())
					rcGetLog() log(RCLogCategory RC_LOG_WARNING, "rcBuildContours: Could not find merge target for bad contour %d.", i)
			} else {
				mcont: RCContour& = cset conts[mergeIdx]
				// Merge by closest points.
				ia, ib: Int
				getClosestIndices(mcont verts, mcont nverts, cont verts, cont nverts, ia, ib)
				if (!mergeContours(mcont, cont, ia, ib)) {
					if (rcGetLog())
						rcGetLog() log(RCLogCategory RC_LOG_WARNING, "rcBuildContours: Failed to merge contours %d and %d.", i, mergeIdx)
				}
			}
		}
	}
	
	simplifyEndTime := rcGetPerformanceTimer()
	endTime := rcGetPerformanceTimer()
	
//	if (rcGetLog())
//	{
//		rcGetLog() log(RC_LOG_PROGRESS, "Create contours: %.3f ms", rcGetDeltaTimeUsec(startTime, endTime)/1000.0f)
//		rcGetLog() log(RC_LOG_PROGRESS, " - boundary: %.3f ms", rcGetDeltaTimeUsec(boundaryStartTime, boundaryEndTime)/1000.0f)
//		rcGetLog() log(RC_LOG_PROGRESS, " - contour: %.3f ms", rcGetDeltaTimeUsec(contourStartTime, contourEndTime)/1000.0f)
//	}
	
	if (rcGetBuildTimes()) {
		rcGetBuildTimes() buildContours += rcGetDeltaTimeUsec(startTime, endTime)
		rcGetBuildTimes() buildContoursTrace += rcGetDeltaTimeUsec(traceStartTime, traceEndTime)
		rcGetBuildTimes() buildContoursSimplify += rcGetDeltaTimeUsec(simplifyStartTime, simplifyEndTime)
	}
	
	return true
}
