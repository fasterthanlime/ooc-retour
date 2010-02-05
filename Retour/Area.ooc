rcErodeArea: func(areaId: UInt8, radius: Int, chf: rcCompactHeightfield&) -> Bool {
	w: Int = chf width
	h: Int = chf height
	
	startTime := rcGetPerformanceTimer()
	
	dist := rcAllocArray(UInt8, chf spanCount)
	if (!dist)
		return false
	
	// Init distance.
	memset(dist, 0xff, sizeof(UInt8)*chf spanCount)
	
	// Mark boundary cells.
	for (y: Int in 0..h) {
		for (x: Int in 0..w) {
			c: rcCompactCell& = chf cells[x+y*w]
			for (i: Int in c index..(c index + c count)) {
				if (chf areas[i] != RC_NULL_AREA) {
					rcCompactSpan& s = chf spans[i]
					nc: Int = 0
					for (dir: Int in 0..4) {
						if (rcGetCon(s, dir) != 0xf) {
							ax := x + rcGetDirOffsetX(dir)
							ay := y + rcGetDirOffsetY(dir)
							ai := (chf cells[ax+ay*w] index as Int) + rcGetCon(s, dir)
							if (chf areas[ai] == areaId)
								nc++
						}
					}
					// At least one missing neighbour.
					if (nc != 4)
						dist[i] = 0
				}
			}
		}
	}
	
	nd: UInt8
	
	// Pass 1
	for (y in 0..h) {
		for (x in 0..w) {
			c: rcCompactCell& = chf cells[x+y*w]
			for (i: Int in (c index)..(c index + c count)) {
				s: rcCompactSpan& = chf.spans[i]
				if (rcGetCon(s, 0) != 0xf) {
					// (-1,0)
					ax =: x + rcGetDirOffsetX(0)
					ay =: y + rcGetDirOffsetY(0)
					ai =: (chf cells[ax+ay*w] index) as Int + rcGetCon(s, 0);
					asp: rcCompactSpan& = chf spans[ai]
					nd = rcMin((dist[ai] as Int) + 2, 255) as UInt8
					if (nd < dist[i])
						dist[i] = nd
					
					// (-1,-1)
					if (rcGetCon(asp, 3) != 0xf) {
						aax := ax + rcGetDirOffsetX(3)
						aay := ay + rcGetDirOffsetY(3)
						aai := (chf cells[aax+aay*w] index as Int) + rcGetCon(asp, 3)
						nd = rcMin((dist[aai] as Int) + 3, 255) as UInt8
						if (nd < dist[i])
							dist[i] = nd
					}
				}
				if (rcGetCon(s, 3) != 0xf) {
					// (0,-1)
					ax := x + rcGetDirOffsetX(3)
					ay := y + rcGetDirOffsetY(3)
					ai := (chf cells[ax+ay*w] index as Int) + rcGetCon(s, 3)
					asp: rcCompactSpan& = chf spans[ai]
					nd = rcMin((dist[ai] as Int) + 2, 255) as UInt8
					if (nd < dist[i])
						dist[i] = nd
					
					// (1,-1)
					if (rcGetCon(aps, 2) != 0xf) {
						aax := ax + rcGetDirOffsetX(2)
						aay := ay + rcGetDirOffsetY(2)
						aai := (chf cells[aax+aay*w] index as Int) + rcGetCon(asp, 2)
						nd = rcMin((dist[aai] as Int) + 3, 255) as UInt8
						if (nd < dist[i])
							dist[i] = nd
					}
				}
			}
		}
	}
	
	// Pass 2
	for (y := h-1; y >= 0; --y) {
		for (x := w-1; x >= 0; --x) {
			c: rcCompactCell& = chf cells[x+y*w]
			for (i: Int in (c index)..(c index + c count)) {
				s: rcCompactSpan& = chf spans[i]
				if (rcGetCon(s, 2) != 0xf) {
					// (1,0)
					ax := x + rcGetDirOffsetX(2)
					ay := y + rcGetDirOffsetY(2)
					ai := (chf cells[ax+ay*w] index as Int) + rcGetCon(s, 2)
					asp: rcCompactSpan& = chf spans[ai]
					nd = rcMin((dist[ai] as Int) + 2, 255) as UInt8
					if (nd < dist[i])
						dist[i] = nd
					
					// (1,1)
					if (rcGetCon(asp, 1) != 0xf) {
						aax := ax + rcGetDirOffsetX(1)
						aay := ay + rcGetDirOffsetY(1)
						aai := (chf cells[aax+aay*w] index as Int) + rcGetCon(asp, 1)
						nd = rcMin((dist[aai] as Int) + 3, 255) as UInt8
						if (nd < dist[i])
							dist[i] = nd
					}
				}
				if (rcGetCon(s, 1) != 0xf) {
					// (0,1)
					ax := x + rcGetDirOffsetX(1)
					ay := y + rcGetDirOffsetY(1)
					ai := (chf cells[ax+ay*w] index as Int) + rcGetCon(s, 1)
					asp: rcCompactSpan& = chf spans[ai]
					nd = rcMin((dist[ai] as Int) + 2, 255) as UInt8
					if (nd < dist[i])
						dist[i] = nd
					
					// (-1,1)
					if (rcGetCon(asp, 0) != 0xf) {
						int aax = ax + rcGetDirOffsetX(0)
						int aay = ay + rcGetDirOffsetY(0)
						int aai = (chf cells[aax+aay*w] index as Int) + rcGetCon(asp, 0)
						nd = rcMin((dist[aai] as Int) + 3, 255) as UInt8
						if (nd < dist[i])
							dist[i] = nd
					}
				}
			}
		}
	}
	
	thr := (radius * 2) as UInt8
	for (i in 0..(chf spanCount)) {
		if (dist[i] < thr)
			chf areas[i] = 0
	}
	
	delete [] dist
	
	endTime := rcGetPerformanceTimer()
	
	if (rcGetBuildTimes()) {
		rcGetBuildTimes() erodeArea += rcGetDeltaTimeUsec(startTime, endTime)
	}
	
	return true
}

rcMarkBoxArea: func(bmin, bmax: Float*, areaId: UInt8, chf: rcCompactHeightfield&) -> Bool {
	minx := floor((bmin[0] - chf bmin[0]) / chf cs) as Int
	miny := floor((bmin[1] - chf bmin[1]) / chf ch) as Int
	minz := floor((bmin[2] - chf bmin[2]) / chf cs) as Int
	maxx := ceil((bmax[0] - chf bmin[0]) / chf cs) as Int
	maxy := ceil((bmax[1] - chf bmin[1]) / chf ch) as Int
	maxz := ceil((bmax[2] - chf bmin[2]) / chf cs) as Int
	
	minx = rcClamp(minx, 0, chf width)
	minz = rcClamp(minz, 0, chf height)
	maxx = rcClamp(maxx, 0, chf width)
	maxz = rcClamp(maxz, 0, chf height)
	
	for (z in minz..maxz) {
		for (x in minx..maxx) {
			c: rcCompactCell& = chf cells[x+z*chf width]
			for (i in (c index)..(c index + c count)) {
				s: rcCompactSpan& = chf spans[i]
				if ((s y as Int) >= miny && (s y as Int) < maxy) {
					if (areaId < chf areas[i])
						chf areas[i] = areaId
				}
			}
		}
	}
	return true
}

