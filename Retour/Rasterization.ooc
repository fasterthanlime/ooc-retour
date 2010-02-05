overlapBounds: inline func(amin, amax, bmin, bmax: Float*) -> Bool {
	overlap := true
	overlap = (amin[0] > bmax[0] || amax[0] < bmin[0]) ? false : overlap
	overlap = (amin[1] > bmax[1] || amax[1] < bmin[1]) ? false : overlap
	overlap = (amin[2] > bmax[2] || amax[2] < bmin[2]) ? false : overlap
	return overlap
}

overlapInterval: inline func(UShort amin, UShort amax, UShort bmin, UShort bmax) -> Bool {
	if (amax < bmin) return false
	if (amin > bmax) return false
	return true
}

allocSpan: static func(hf: RCHeightfield&) -> RCSpan* {
	// If running out of memory, allocate new page and update the freelist.
	if (!hf freelist || !hf freelist next) {
		// Create new page.
		// Allocate memory for the new pool.
		size := (sizeof(RCSpanPool)-sizeof(RCSpan)) + sizeof(RCSpan)*RC_SPANS_PER_POOL
		pool := rcAllocArray(UInt8, size) as RCSpanPool*
		if (!pool) return 0
		pool next = 0
		// Add the pool Into the list of pools.
		pool next = hf pools
		hf pools = pool
		// Add new items to the free list.
		freelist := hf freelist
		head: RCSpan* = &pool items[0]
		it: RCSpan* = &pool items[RC_SPANS_PER_POOL]
		while (it != head) { // do
			--it
			it next = freelist
			freelist = it
		} //while (it != head)
		hf freelist = it
	}
	
	// Pop item from in front of the free list.
	it: RCSpan* = hf freelist
	hf freelist = hf freelist next
	return it
}

freeSpan: static func(RCHeightfield& hf, RCSpan* ptr) {
	if (!ptr) return
	// Add the node in front of the free list.
	ptr next = hf freelist
	hf freelist = ptr
}

addSpan: static func(hf: RCHeightfield&, x, y: Int, smin, smax, flags: UShort, flagMergeThr: Int) {
	idx := x + y*hf width
	
	s := allocSpan(hf)
	s smin = smin
	s smax = smax
	s flags = flags
	s next = 0
	
	// Empty cell, add he first span.
	if (!hf spans[idx]) {
		hf spans[idx] = s
		return
	}
	prev: RCSpan* = 0
	cur := hf spans[idx]
	
	// Insert and merge spans.
	while (cur) {
		if (cur smin > s smax) {
			// Current span is further than the new span, break.
			break
		} else if (cur smax < s smin) {
			// Current span is before the new span advance.
			prev = cur
			cur = cur next
		} else {
			// Merge spans.
			if (cur smin < s smin)
				s smin = cur smin
			if (cur smax > s smax)
				s smax = cur smax
			
			// Merge flags.
			if (rcAbs((s smax as Int) - (cur smax as Int)) <= flagMergeThr)
				s flags |= cur flags
			
			// Remove current span.
			next := cur next
			freeSpan(hf, cur)
			if (prev)
				prev next = next
			else
				hf spans[idx] = next
			cur = next
		}
	}
	
	// Insert new span.
	if (prev) {
		s next = prev next
		prev next = s
	} else {
		s next = hf spans[idx]
		hf spans[idx] = s
	}
}

clipPoly: static func(inz: Float*, n: Int, out: Float*, pnx, pnz, pd: Float) -> Int {
	d: Float[12]
	for (i: Int in 0..n)
		d[i] = pnx*inz[i*3+0] + pnz*inz[i*3+2] + pd
	
	m := 0
	j := n-1
	for (i: Int in 0..n) {
		ina := d[j] >= 0
		inb := d[i] >= 0
		if (ina != inb) {
			s: Float = d[j] / (d[j] - d[i])
			out[m*3+0] = in[j*3+0] + (in[i*3+0] - in[j*3+0])*s
			out[m*3+1] = in[j*3+1] + (in[i*3+1] - in[j*3+1])*s
			out[m*3+2] = in[j*3+2] + (in[i*3+2] - in[j*3+2])*s
			m++
		}
		if (inb) {
			out[m*3+0] = in[i*3+0]
			out[m*3+1] = in[i*3+1]
			out[m*3+2] = in[i*3+2]
			m++
		}
		j=i
	}
	return m
}

rasterizeTri: static func(v0, v1, v2: Float*, flags: UInt8, hf: RCHeightfield&, bmin, bmax: Float*, cs, ics, ich: Float, flagMergeThr: Int) {
	w := hf width
	h := hf height
	tmin, tmax: Float[3]
	by: Float = bmax[1] - bmin[1]
	
	// Calculate the bounding box of the triangle.
	vcopy(tmin, v0)
	vcopy(tmax, v0)
	vmin(tmin, v1)
	vmin(tmin, v2)
	vmax(tmax, v1)
	vmax(tmax, v2)
	
	// If the triangle does not touch the bbox of the heightfield, skip the triagle.
	if (!overlapBounds(bmin, bmax, tmin, tmax))
		return
	
	// Calculate the footpring of the triangle on the grid.
	x0 := ((tmin[0] - bmin[0])*ics) as Int
	y0 := ((tmin[2] - bmin[2])*ics) as Int
	x1 := ((tmax[0] - bmin[0])*ics) as Int
	y1 := ((tmax[2] - bmin[2])*ics) as Int
	x0 = rcClamp(x0, 0, w-1)
	y0 = rcClamp(y0, 0, h-1)
	x1 = rcClamp(x1, 0, w-1)
	y1 = rcClamp(y1, 0, h-1)
	
	// Clip the triangle Into all grid cells it touches.
	inz, out, inrow: Float[7*3]
	
	for (y := y0; y <= y1; ++y) {
		// Clip polygon to row.
		vcopy(&inz[0], v0)
		vcopy(&inz[1*3], v1)
		vcopy(&inz[2*3], v2)
		nvrow := 3
		cz: Float = bmin[2] + y*cs
		nvrow = clipPoly(inz, nvrow, out, 0, 1, -cz)
		if (nvrow < 3) continue
		nvrow = clipPoly(out, nvrow, inrow, 0, -1, cz+cs)
		if (nvrow < 3) continue
		
		for (x := x0; x <= x1; ++x) {
			// Clip polygon to column.
			nv := nvrow
			cx := bmin[0] + x*cs
			nv = clipPoly(inrow, nv, out, 1, 0, -cx)
			if (nv < 3) continue
			nv = clipPoly(out, nv, inz, -1, 0, cx+cs)
			if (nv < 3) continue
			
			// Calculate min and max of the span.
			smin: Float = inz[1], smax = inz[1]
			for (i: Int in 1..nv) {
				smin = rcMin(smin, inz[i*3+1])
				smax = rcMax(smax, inz[i*3+1])
			}
			smin -= bmin[1]
			smax -= bmin[1]
			// Skip the span if it is outside the heightfield bbox
			if (smax < 0.0) continue
			if (smin > by) continue
			// Clamp the span to the heightfield bbox.
			if (smin < 0.0) smin = 0
			if (smax > by) smax = by
			
			// Snap the span to the heightfield height grid.
			ismin := rcClamp(floor(smin * ich) as Int, 0, 0x7fff) as UShort
			ismax := rcClamp(ceil(smax * ich) as Int, 0, 0x7fff) as UShort
			
			addSpan(hf, x, y, ismin, ismax, flags, flagMergeThr)
		}
	}
}

rcRasterizeTriangle: func(v0, v1, v2: Float*, flags: UInt8, solid: RCHeightfield&, flagMergeThr: Int) {
	startTime := rcGetPerformanceTimer()

	ics: Float = 1.0/solid cs
	ich: Float = 1.0/solid ch
	rasterizeTri(v0, v1, v2, flags, solid, solid bmin, solid bmax, solid cs, ics, ich, flagMergeThr)
	endTime := rcGetPerformanceTimer()
	if (rcGetBuildTimes())
		rcGetBuildTimes() rasterizeTriangles += rcGetDeltaTimeUsec(startTime, endTime)
}

rcRasterizeTriangles: func(verts: Float*, nv: Int, tris: Int*, flags: UInt8, nt: Int, solid: RCHeightfield&, flagMergeThr: Int) {
	startTime := rcGetPerformanceTimer()
	
	ics = 1.0/solid cs
	ich = 1.0/solid ch
	// Rasterize triangles.
	for (i: Int in 0..nt) {
		v0: Float* = &verts[tris[i*3+0]*3]
		v1: Float* = &verts[tris[i*3+1]*3]
		v2: Float* = &verts[tris[i*3+2]*3]
		// Rasterize.
		rasterizeTri(v0, v1, v2, flags[i], solid, solid bmin, solid bmax, solid cs, ics, ich, flagMergeThr)
	}
	endTime := rcGetPerformanceTimer()
	if (rcGetBuildTimes())
		rcGetBuildTimes() rasterizeTriangles += rcGetDeltaTimeUsec(startTime, endTime)
}

rcRasterizeTriangles: func(verts: Float*, nv: Int, tris: UShort*, flags: UInt8*, nt: Int, solid: RCHeightfield&, flagMergeThr: Int) {
	startTime := rcGetPerformanceTimer()
	
	ics: Float = 1.0/solid cs
	ich: Float = 1.0/solid ch
	// Rasterize triangles.
	for (i: Int in 0..nt) {
		v0: Float* = &verts[tris[i*3+0]*3]
		v1: Float* = &verts[tris[i*3+1]*3]
		v2: Float* = &verts[tris[i*3+2]*3]
		// Rasterize.
		rasterizeTri(v0, v1, v2, flags[i], solid, solid bmin, solid bmax, solid cs, ics, ich, flagMergeThr)
	}
	endTime := rcGetPerformanceTimer()
	if (rcGetBuildTimes())
		rcGetBuildTimes() rasterizeTriangles += rcGetDeltaTimeUsec(startTime, endTime)
}

rcRasterizeTriangles func(verts: Float*, flags: UInt8*, nt: Int, solid: RCHeightfield&, flagMergeThr: Int) {
	startTime := rcGetPerformanceTimer()
	
	ics: Float = 1.0/solid cs
	ich: Float = 1.0/solid ch
	// Rasterize triangles.
	for (i: Int in 0..nt) {
		v0 = &verts[(i*3+0)*3]
		v1 = &verts[(i*3+1)*3]
		v2 = &verts[(i*3+2)*3]
		// Rasterize.
		rasterizeTri(v0, v1, v2, flags[i], solid, solid bmin, solid bmax, solid cs, ics, ich, flagMergeThr)
	}
	endTime := rcGetPerformanceTimer()
	if (rcGetBuildTimes())
		rcGetBuildTimes() rasterizeTriangles += rcGetDeltaTimeUsec(startTime, endTime)
}

