use retour
import structs/Array
import lang/math
import Retour/[Common, Log, Timer] //, Area, Contour, Filter, Mesh, MeshDetail, Rasterization, Region

RCConfig: class {
	init: func() {}
	
	width, height: Int				// Dimensions of the rasterized heighfield (vx)
	tileSize: Int					// Width and Height of a tile (vx)
	borderSize: Int					// Non-navigable Border around the heightfield (vx)
	cs, ch: Float					// Grid cell size and height (wu)
	bmin, bmax: Float[3]			// Grid bounds (wu)
	walkableSlopeAngle: Float		// Maximum walkble slope angle in degrees.
	walkableHeight: Int				// Minimum height where the agent can still walk (vx)
	walkableClimb: Int				// Maximum height between grid cells the agent can climb (vx)
	walkableRadius: Int				// Radius of the agent in cells (vx)
	maxEdgeLen: Int					// Maximum contour edge length (vx)
	maxSimplificationError: Float	// Maximum distance error from contour to cells (vx)
	minRegionSize: Int				// Minimum regions size. Smaller regions will be deleted (vx)
	mergeRegionSize: Int			// Minimum regions size. Smaller regions will be merged (vx)
	maxVertsPerPoly: Int			// Max number of vertices per polygon
	detailSampleDist: Float			// Detail mesh sample spacing.
	detailSampleMaxError: Float		// Detail mesh simplification max sample error.
}

RCSpan: class {
	init: func() {}
	
	smin: UInt //:15		// Span min height.
	smax: UInt //:15		// Span max height.
	flags: UInt //:2		// Span flags.
	next: RCSpan*			// Next span in column.
}

RC_SPANS_PER_POOL: static const Int = 2048

RCSpanPool: class {
	init: func() {}
	
	next: RCSpanPool*	// Pointer to next pool.
	items: RCSpan[1]	// Array of spans (size RC_SPANS_PER_POOL).
}

RCHeightfield: class {
	init: func() {}
	
	width, height: Int			// Dimension of the heightfield.
	bmin, bmax: Float[3]		// Bounding box of the heightfield
	cs, ch: Float				// Cell size and height.
	spans: Array<RCSpan>		// Heightfield of spans (width*height).
	pools: RCSpanPool*			// Linked list of span pools.
	freelist: RCSpan*			// Pointer to next free span.
}

RCCompactCell: class {
	init: func() {}
	index: UInt //:24	// Index to first span in column.
	count: UInt //:8	// Number of spans in this column.
}

RCCompactSpan: class {
	init: func() {}
	y: UShort		// Bottom coordinate of the span.
	con: UShort		// Connections to neighbour cells.
	h: UInt8		// Height of the span.
}

RCCompactHeightfield: class {
	init: func() {}
	
	width, height: Int					// Width and height of the heighfield.
	spanCount: Int						// Number of spans in the heightfield.
	walkableHeight, walkableClimb: Int	// Agent properties.
	maxDistance: UShort			// Maximum distance value stored in heightfield.
	maxRegions: UShort			// Maximum Region Id stored in heightfield.
	bmin, bmax: Float[3]				// Bounding box of the heightfield.
	cs, ch: Float						// Cell size and height.
	cells: Array<RCCompactCell>				// Pointer to width*height cells.
	spans: Array<RCCompactSpan>				// Pointer to spans.
	dist: Array<UShort>				// Pointer to per span distance to border.
	regs: Array<UShort>				// Pointer to per span region ID.
	areas: Array<UInt8>				// Pointer to per span area ID.
}

RCContour: class {
	init: func() {}
	
	verts: Array<Int>			// Vertex coordinates, each vertex contains 4 components.
	nverts: Int			// Number of vertices.
	rverts: Array<Int>		// Raw vertex coordinates, each vertex contains 4 components.
	nrverts: Int		// Number of raw vertices.
	reg: UShort			// Region ID of the contour.
	area: UInt8			// Area ID of the contour.
}

RCContourSet: class {
	init: func() {}
	
	conts: RCContour* 		// Pointer to all contours.
	nconts: Int				// Number of contours.
	bmin, bmax: Float[3]	// Bounding box of the heightfield.
	cs, ch: Float			// Cell size and height.
}

// Polymesh stores a connected mesh of polygons.
// The polygons are store in an array where each polygons takes
// 'nvp*2' elements. The first 'nvp' elements are indices to vertices
// and the second 'nvp' elements are indices to neighbour polygons.
// If a polygona has less than 'bvp' vertices, the remaining indices
// are set to RC_MESH_NULL_IDX. If an polygon edge does not have a neighbour
// the neighbour index is set to RC_MESH_NULL_IDX.
// Vertices can be transformed into world space as follows:
//   x = bmin[0] + verts[i*3+0]*cs
//   y = bmin[1] + verts[i*3+1]*ch
//   z = bmin[2] + verts[i*3+2]*cs
RCPolyMesh: class {
	init: func() {
		nvp = 3
	}
	
	verts: Array<UShort>			// Vertices of the mesh, 3 elements per vertex.
	polys: Array<UShort>			// Polygons of the mesh, nvp*2 elements per polygon.
	regs: Array<UShort>			// Region ID of the polygons.
	areas: Array<UInt8>			// Area ID of polygons.
	nverts: Int			// Number of vertices.
	npolys: Int			// Number of polygons.
	nvp: Int			// Max number of vertices per polygon.
	bmin, bmax: Float[3]	// Bounding box of the mesh.
	cs, ch: Float			// Cell size and height.
}

// Detail mesh generated from a RCPolyMesh.
// Each submesh represents a polygon in the polymesh and they are stored in
// excatly same order. Each submesh is described as 4 values:
// base vertex, vertex count, base triangle, triangle count. That is,
//   UInt8* t = &dtl.tris[(tbase+i)*3]; and
//   float* v = &dtl.verts[(vbase+t[j])*3]
// If the input polygon has 'n' vertices, those vertices are first in the
// submesh vertex list. This allows to compres the mesh by not storing the
// first vertices and using the polymesh vertices instead.
RCPolyMeshDetail: class {
	init: func() {}
	
	meshes: Array<UShort>		// Pointer to all mesh data.
	verts: Array<Float>		// Pointer to all vertex data.
	tris: Array<UInt8>		// Pointer to all triangle data.
	nmeshes: Int		// Number of meshes.
	nverts: Int			// Number of total vertices.
	ntris: Int			// Number of triangles.
}

RCIntArray: class {
	m_data: Int[]
	m_size, m_cap: Int
	
	init: func() {}
	init: func~withSize(size: Int) {
		m_data = rcAllocArray(Int, size)
		m_cap = size
	}
	
	resize: func(n: Int) {
		if (n > m_cap) {
			if (!m_cap) m_cap = 8
			while (m_cap < n) m_cap *= 2
			newData := rcAllocArray(Int, m_cap)
			if (m_size && newData) memcpy(newData, m_data, m_size*sizeof(Int))
			//delete [] m_data
			m_data = newData
		}
		m_size = n
	}
	
	push: func(item: Int) {
		resize(m_size + 1)
		m_data[m_size - 1] = item
	}
	
	pop: func() -> Int {
		if (m_size > 0)
			m_size -= 1
		return m_data[m_size]
	}
	
	size: func() -> Int {
		return m_size
	}
}

operator [] (array: RCIntArray, index: Int) -> Int {
	return array m_data[index]
}

operator []= (array: RCIntArray, index: Int, value: Int) {
	// Might need to use push.. there was no []= operator defined
	// in the C++ code, yet it uses iarray[blah] = blah
	array m_data[index] = value
}

/*
RCScopedDelete: class <T> {
	ptr: T*
	
	init: func() {}
	init: func~withPtr(=ptr)
	destroy: func() {
		delete [] ptr
	}
	
	inline operator T*() { return ptr; }
	inline T* operator=(T* p) { ptr = p; return ptr; }
}*/

// enum
RCSpanFlags: class {
	RC_WALKABLE: static const Int = 0x01
	RC_LEDGE: static const Int = 0x02
}

// If heightfield region ID has the following bit set, the region is on border area
// and excluded from many calculations.
RC_BORDER_REG: static const UShort = 0x8000

// If contour region ID has the following bit set, the vertex will be later
// removed in order to match the segments and vertices at tile boundaries.
RC_BORDER_VERTEX: static const  Int = 0x10000

RC_AREA_BORDER: static const  Int = 0x20000

// Mask used with contours to extract region id.
RC_CONTOUR_REG_MASK: static const  Int = 0xffff

// Null index which is used with meshes to mark unset or invalid indices.
RC_MESH_NULL_IDX: static const  UShort = 0xffff

// Area ID that is considered empty.
RC_NULL_AREA: static const  UInt8 = 0

// Area ID that is considered generally walkable.
RC_WALKABLE_AREA: static const  UInt8 = 255

// Value returned by rcGetCon() if the direction is not connected.
RC_NOT_CONNECTED: static const  Int = 0xf

// Compact span neighbour helpers.
rcSetCon: inline func(s: RCCompactSpan@, dir, i: Int) {
	s con &= ~(0xf << (dir*4))
	s con |= (i & 0xf) << (dir*4)
}

rcGetCon: inline func(s: RCCompactSpan@, dir: Int) -> Int {
	return (s con >> (dir*4)) & 0xf
}

rcGetDirOffsetX: inline func(dir: Int) -> Int {
	offset: Int[4] = [-1, 0, 1, 0]
	return offset[dir & 0x03]
}

rcGetDirOffsetY: inline func(dir: Int) -> Int {
	offset: Int[4] = [0, 1, 0, -1]
	return offset[dir & 0x03]
}

rcCalcBounds: func(verts: Float*, nv: Int, bmin, bmax: Float*) {
	// Calculate bounding box.
	vcopy(bmin, verts)
	vcopy(bmax, verts)
	for (i in 1..nv) {
		v: Float* = verts[i*3]&
		vmin(bmin, v)
		vmax(bmax, v)
	}
}

rcCalcGridSize: func(bmin, bmax: Float*, cs: Float, w, h: Int*) {
	w@ = ((bmax[0] - bmin[0])/cs+0.5) as Int
	h@ = ((bmax[2] - bmin[2])/cs+0.5) as Int
}

rcCreateHeightfield: func(hf: RCHeightfield@, width, height: Int, bmin, bmax: Float*, cs, ch: Float) -> Bool {
	hf width = width
	hf height = height
	hf spans = Array<RCSpan> new(hf width * hf height)
	vcopy(hf bmin, bmin)
	vcopy(hf bmax, bmax)
	hf cs = cs
	hf ch = ch
	if (!hf spans)
		return false
	//memset(hf spans data, 0, hf spans size)
	return true
}

calcTriNormal: static func(v0, v1, v2: Float*, norm: Float*) {
	e0, e1: Float[3]
	vsub(e0, v1, v0)
	vsub(e1, v2, v0)
	vcross(norm, e0, e1)
	vnormalize(norm)
}

rcMarkWalkableTriangles: func(walkableSlopeAngle: Float, verts: Float*, nv: Int, tris: Int*, nt: Int, flags: UInt8*) {
	walkableThr: Float = cos(walkableSlopeAngle/180.0 * RCConstants PI)
	norm: Float[3]
	
	for (i in 0..nt) {
		tri: Int* = tris[i*3]&
		calcTriNormal(verts[tri[0]*3]&, verts[tri[1]*3]&, verts[tri[2]*3]&, norm)
		// Check if the face is walkable.
		if (norm[1] > walkableThr)
			flags[i] |= RCSpanFlags RC_WALKABLE
	}
}

getSpanCount: static func(flags: UInt8, hf: RCHeightfield@) -> Int {
	w := hf width
	h := hf height
	spanCount := 0
	for (y: Int in 0..h) {
		for (x: Int in 0..w) {
			s := hf spans[x + y*w]
			while (s) {
				if (s flags == flags)
					spanCount += 1
				s = s next
			}
		}
	}
	return spanCount
}

rcBuildCompactHeightfield: func(walkableHeight, walkableClimb: Int, flags: UInt8, hf: RCHeightfield@, chf: RCCompactHeightfield@) -> Bool {
	startTime := rcGetPerformanceTimer()
	
	w := hf width
	h := hf height
	spanCount := getSpanCount(flags, hf)
	
	// Fill in header.
	chf width = w
	chf height = h
	chf spanCount = spanCount
	chf walkableHeight = walkableHeight
	chf walkableClimb = walkableClimb
	chf maxRegions = 0
	vcopy(chf bmin, hf bmin)
	vcopy(chf bmax, hf bmax)
	chf bmax[1] += walkableHeight * hf ch
	chf cs = hf cs
	chf ch = hf ch
	
	chf cells = Array<RCCompactCell> new(w*h)
	if (!chf cells) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildCompactHeightfield: Out of memory 'chf cells' (%d)", w*h)
		return false
	}
	//memset(chf cells, 0, sizeof(RCCompactCell)*w*h)
	
	chf spans = Array<RCCompactSpan> new(spanCount)
	if (!chf spans) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildCompactHeightfield: Out of memory 'chf spans' (%d)", spanCount)
		return false
	}
	//memset(chf spans, 0, sizeof(RCCompactSpan)*spanCount)
	
	chf areas = Array<UInt8> new(spanCount)
	if (!chf areas) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildCompactHeightfield: Out of memory 'chf areas' (%d)", spanCount)
		return false
	}
	memset(chf areas data, RC_WALKABLE_AREA, chf areas size)
	
	MAX_HEIGHT: static const Int = 0xffff
	
	// Fill in cells and spans.
	idx := 0
	for (y: Int in 0..h) {
		for (x: Int in 0..w) {
			s := hf spans[x + y*w]
			// If there are no spans at this cell, just leave the data to index=0, count=0.
			if (!s) continue
			c := chf cells[x+y*w]
			c index = idx
			c count = 0
			while (s) {
				if (s flags == flags) {
					bot := (s smax as Int)
					top := s next ? (s next smin as Int) : MAX_HEIGHT
					chf spans[idx] y = rcClamp(bot, 0, 0xffff) as UShort
					chf spans[idx] h = rcClamp(top - bot, 0, 0xff) as UInt8
					idx += 1
					c count += 1
				}
				s = s next
			}
		}
	}
	
	// Find neighbour connections.
	for (y: Int in 0..h) {
		for (x: Int in 0..w) {
			 c := chf cells[x+y*w]
			for (i: Int in (c index)..(c index + c count)) {
				s := chf spans[i]
				for (dir: Int in 0..4) {
					rcSetCon(s, dir, RC_NOT_CONNECTED)
					nx: Int = x + rcGetDirOffsetX(dir)
					ny: Int = y + rcGetDirOffsetY(dir)
					// First check that the neighbour cell is in bounds.
					if (nx < 0 || ny < 0 || nx >= w || ny >= h)
						continue
						
					// Iterate over all neighbour spans and check if any of the is
					// accessible from current cell.
					nc := chf cells[nx+ny*w]
					for (k: Int in (nc index)..(nc index+nc count)) {
						ns := chf spans[k]
						bot := rcMax(s y, ns y)
						top := rcMin(s y+s h, ns y+ns h)

						// Check that the gap between the spans is walkable,
						// and that the climb height between the gaps is not too high.
						if ((top - bot) >= walkableHeight && rcAbs((ns y as Int) - (s y as Int)) <= walkableClimb) {
							// Mark direction as walkable.
							rcSetCon(s, dir, k - (nc index as Int))
							break
						}
					}
				}
			}
		}
	}
	
	endTime := rcGetPerformanceTimer()
	
	if (rcGetBuildTimes())
		rcGetBuildTimes() buildCompact += rcGetDeltaTimeUsec(startTime, endTime)
	
	return true
}

getHeightfieldMemoryUsage: static func(hf: RCHeightfield@) -> Int {
	size := 0
	size += sizeof(hf)
	size += hf width * hf height * sizeof(Pointer) // sizeof(RCSpan*)
	
	pool: RCSpanPool* = hf pools
	while (pool) {
		size += (sizeof(RCSpanPool) - sizeof(RCSpan)) + sizeof(RCSpan)*RC_SPANS_PER_POOL
		pool = pool next
	}
	return size
}

getCompactHeightFieldMemoryusage: static func(chf: RCCompactHeightfield@) -> Int {
	size := 0
	size += sizeof(RCCompactHeightfield)
	size += sizeof(RCCompactSpan) * chf spanCount
	size += sizeof(RCCompactCell) * chf width * chf height
	return size
}

