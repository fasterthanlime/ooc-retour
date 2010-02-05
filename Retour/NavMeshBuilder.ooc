
DTNavMeshCreateParams: class {
	// Navmesh vertices.
	verts: UShort*
	vertCount: Int
	// Navmesh polygons
	polys: UShort*
	polyCount: Int
	nvp: Int
	// Navmesh Detail
	detailMeshes: UShort*
	detailVerts: Float*
	detailVertsCount: Int
	detailTris: UInt8*
	detailTriCount: Int
	// Off-Mesh Connections.
	offMeshConVerts: Float*
	offMeshConRad: Float*
	offMeshConDir: UInt8*
	offMeshConCount: Int
	// Settings
	walkableHeight: Float
	walkableRadius: Float
	walkableClimb: Float
	bmin, bmax: Float[3]
	cs: Float
	ch: Float
	tileSize: Int
}

MESH_NULL_IDX: static const UShort = 0xffff

BVItem: class {
	bmin, bmax: UShort[3]
	i: Int
}

compareItemX: static func(va, vb: Void*) -> Int {
	a := (va as BVItem*)
	b := (vb as BVItem*)
	if (a bmin[0] < b bmin[0])
		return -1
	if (a bmin[0] > b bmin[0])
		return 1
	return 0
}

compareItemY: static func(va, vb: Void*) -> Int {
	a := (va as BVItem*)
	b := (vb as BVItem*)
	if (a bmin[1] < b bmin[1])
		return -1
	if (a bmin[1] > b bmin[1])
		return 1
	return 0
}

compareItemZ: static func(va, vb: Void*) -> Int {
	a := (va as BVItem*)
	b := (vb as BVItem*)
	if (a bmin[2] < b bmin[2])
		return -1
	if (a bmin[2] > b bmin[2])
		return 1
	return 0
}

calcExtends: static func(items: BVItem*, nitems, imin, imax: Int, bmin, bmax: UShort*) {
	bmin[0] = items[imin] bmin[0]
	bmin[1] = items[imin] bmin[1]
	bmin[2] = items[imin] bmin[2]
	
	bmax[0] = items[imin] bmax[0]
	bmax[1] = items[imin] bmax[1]
	bmax[2] = items[imin] bmax[2]
	
	for (i: Int in (imin+1)..imax) {
		it: BVItem& = items[i]
		if (it bmin[0] < bmin[0]) bmin[0] = it bmin[0]
		if (it bmin[1] < bmin[1]) bmin[1] = it bmin[1]
		if (it bmin[2] < bmin[2]) bmin[2] = it bmin[2]
		
		if (it bmax[0] > bmax[0]) bmax[0] = it bmax[0]
		if (it bmax[1] > bmax[1]) bmax[1] = it bmax[1]
		if (it bmax[2] > bmax[2]) bmax[2] = it bmax[2]
	}
}

longestAxis: inline func(x, y, z: UShort) -> Int {
	axis := 0
	maxVal := x
	if (y > maxVal) {
		axis = 1
		maxVal = y
	}
	if (z > maxVal) {
		axis = 2
		maxVal = z
	}
	return axis
}

subdivide: static func(items: BVItem*, nitems, imin, imax, curNode: Int, nodes: DTBVNode*) {
	inum := imax - imin
	icur: Int = curNode
	
	node: DTBVNode& = nodes[curNode++]
	if (inum == 1) {
		// Leaf
		node bmin[0] = items[imin] bmin[0]
		node bmin[1] = items[imin] bmin[1]
		node bmin[2] = items[imin] bmin[2]
		
		node bmax[0] = items[imin] bmax[0]
		node bmax[1] = items[imin] bmax[1]
		node bmax[2] = items[imin] bmax[2]
		
		node i = items[imin] i
	} else {
		// Split
		calcExtends(items, nitems, imin, imax, node bmin, node bmax)
		axis := longestAxis(node bmax[0] - node bmin[0], node bmax[1] - node bmin[1], node bmax[2] - node bmin[2])
		if (axis == 0) {
			// Sort along x-axis
			qsort(items+imin, inum, sizeof(BVItem), compareItemX)
		} else if (axis == 1) {
			// Sort along y-axis
			qsort(items+imin, inum, sizeof(BVItem), compareItemY)
		} else {
			// Sort along z-axis
			qsort(items+imin, inum, sizeof(BVItem), compareItemZ)
		}
		isplit := imin+inum/2
		
		// Left
		subdivide(items, nitems, imin, isplit, curNode, nodes)
		// Right
		subdivide(items, nitems, isplit, imax, curNode, nodes)
		
		iescape := curNode - icur
		// Negative index means escape.
		node i = -iescape
	}
}

createBVTree: static func(verts: UShort*, nverts: Int, polys: UShort*, npolys, nvp: Int,
						cs, ch: Float, nnodes: Int, nodes: DTBVNode*) -> Int {
	// Build tree
	items := rcAllocArray(BVItem, npolys)
	for (i: Int in 0..npolys) {
		it: BVItem& = items[i]
		it i = i
		// Calc polygon bounds.
		p: UShort* = &polys[i*nvp*2]
		it bmin[0] = it bmax[0] = verts[p[0]*3+0]
		it bmin[1] = it bmax[1] = verts[p[0]*3+1]
		it bmin[2] = it bmax[2] = verts[p[0]*3+2]
		
		for (j: Int in 1..nvp) {
			if (p[j] == MESH_NULL_IDX) break
			x := verts[p[j]*3+0]
			y := verts[p[j]*3+1]
			z := verts[p[j]*3+2]
			
			if (x < it bmin[0]) it bmin[0] = x
			if (y < it bmin[1]) it bmin[1] = y
			if (z < it bmin[2]) it bmin[2] = z
			
			if (x > it bmax[0]) it bmax[0] = x
			if (y > it bmax[1]) it bmax[1] = y
			if (z > it bmax[2]) it bmax[2] = z
		}
		// Remap y
		it bmin[1] = floor((it bmin[1] as Float)*ch/cs) as UShort
		it bmax[1] = ceil((it bmax[1] as Float)*ch/cs) as UShort
	}
	
	curNode := 0
	subdivide(items, npolys, 0, npolys, curNode, nodes)
	delete [] items
	
	return curNode
}

classifyOffMeshPoint: static func(pt, bmin, bmax: Float*) -> UInt8 {
	XP: static const UInt8 = 1<<0
	ZP: static const UInt8 = 1<<1
	XM: static const UInt8 = 1<<2
	ZM: static const UInt8 = 1<<3
	
	outcode: UInt8 = 0
	outcode |= (pt[0] >= bmax[0]) ? XP : 0
	outcode |= (pt[2] >= bmax[2]) ? ZP : 0
	outcode |= (pt[0] < bmin[0])  ? XM : 0
	outcode |= (pt[2] < bmin[2])  ? ZM : 0
	
	match (outcode) {
		case XP: return 0
		case XP|ZP: return 1
		case ZP: return 2
		case XM|ZP: return 3
		case XM: return 4
		case XM|ZM: return 5
		case ZM: return 6
		case XP|ZM: return 7
	}
	return 0xff
}

// TODO: Better error handling.
dtCreateNavMeshData: func(params: DTNavMeshCreateParams*, outData: UInt8**, outDataSize: Int*) -> Bool {
	if (params nvp > DT_VERTS_PER_POLYGON)
		return false
	if (params vertCount >= 0xffff)
		return false
	if (!params vertCount || !params verts)
		return false
	if (!params polyCount || !params polys)
		return false
	if (!params detailMeshes || !params detailVerts || !params detailTris)
		return false
	
	nvp := params nvp
	
	// Classify off-mesh connection points. We store only the connections
	// whose start point is inside the tile.
	offMeshConFlags: UInt8* = new UInt8 [params offMeshConCount*2]
	if (!offMeshConFlags)
		return false
	
	storedOffMeshConCount := 0
	offMeshConLinkCount := 0

	for (i: Int in 0..(params offMeshConCount)) {
		offMeshConFlags[i*2+0] = classifyOffMeshPoint(&params offMeshConVerts[(i*2+0)*3], params bmin, params bmax)
		offMeshConFlags[i*2+1] = classifyOffMeshPoint(&params offMeshConVerts[(i*2+1)*3], params bmin, params bmax)

		// Cound how many links should be allocated for off-mesh connections.
		if (offMeshConFlags[i*2+0] == 0xff)
			offMeshConLinkCount++
		if (offMeshConFlags[i*2+1] == 0xff)
			offMeshConLinkCount++

		if (offMeshConFlags[i*2+0] == 0xff)
			storedOffMeshConCount++
	}
	
	// Off-mesh connectionss are stored as polygons, adjust values.
	totPolyCount := params polyCount + storedOffMeshConCount
	totVertCount := params vertCount + storedOffMeshConCount*2
	
	// Find portal edges which are at tile borders.
	edgeCount := 0
	portalCount := 0
	for (i: Int in 0..(params polyCount)) {
		p: UShort* = &params polys[i*2*nvp]
		for (j: Int in 0..nvp) {
			if (p[j] == MESH_NULL_IDX) break
			nj := j+1
			if (nj >= nvp || p[nj] == MESH_NULL_IDX) nj = 0
			va: UShort* = &params verts[p[j]*3]
			vb: UShort* = &params verts[p[nj]*3]
			edgeCount++
			if (params tileSize > 0) {
				if (va[0] == params tileSize && vb[0] == params tileSize)
					portalCount++; // x+
				else if (va[2] == params tileSize && vb[2] == params tileSize)
					portalCount++; // z+
				else if (va[0] == 0 && vb[0] == 0)
					portalCount++; // x-
				else if (va[2] == 0 && vb[2] == 0)
					portalCount++; // z-
			}
		}
	}
	
	maxLinkCount := edgeCount + portalCount*2 + offMeshConLinkCount*2
	
	// Find unique detail vertices.
	uniqueDetailVertCount := 0
	for (i: Int in 0..(params polyCount)) {
		p: UShort* = &params polys[i*nvp*2]
		ndv := params detailMeshes[i*4+1]
		nv := 0
		for (j: Int in 0..nvp) {
			if (p[j] == MESH_NULL_IDX) break
			nv++
		}
		ndv -= nv
		uniqueDetailVertCount += ndv
	}
	
	// Calculate data size
	headerSize := align4(sizeof(DTMeshHeader))
	vertsSize := align4(sizeof(Float)*3*totVertCount)
	polysSize := align4(sizeof(DTPoly)*totPolyCount)
	linksSize := align4(sizeof(DTLink)*maxLinkCount)
	detailMeshesSize := align4(sizeof(DTPolyDetail)*params polyCount)
	detailVertsSize := align4(sizeof(Float)*3*uniqueDetailVertCount)
	detailTrisSize := align4(sizeof(UInt8)*4*params detailTriCount)
	bvTreeSize := align4(sizeof(DTBVNode)*params polyCount*2)
	offMeshConsSize := align4(sizeof(DTOffMeshConnection)*storedOffMeshConCount)
	
	dataSize := headerSize + vertsSize + polysSize + linksSize + detailMeshesSize + detailVertsSize + detailTrisSize + bvTreeSize + offMeshConsSize
	
	data := rcAllocArray(UInt8, dataSize)
	if (!data)
		return false
	memset(data, 0, dataSize)
	
	d: UInt8* = data
	header: DTMeshHeader* = (d as DTMeshHeader*); d += headerSize
	navVerts: Float* = (d as Float*); d += vertsSize
	navPolys: DTPoly* = (d as DTPoly*); d += polysSize
	d += linksSize
	navDMeshes: DTPolyDetail* = (d as DTPolyDetail*); d += detailMeshesSize
	navDVerts: Float* = (d as Float*); d += detailVertsSize
	navDTris: UInt8* = (d as UInt8*); d += detailTrisSize
	navBvtree: DTBVNode* = (d as DTBVNode*); d += bvTreeSize
	offMeshCons: DTOffMeshConnection* = (DTOffMeshConnection*)d; d += offMeshConsSize
	
	// Store header
	header magic = DT_NAVMESH_MAGIC
	header dversion = DT_NAVMESH_VERSION
	header polyCount = totPolyCount
	header vertCount = totVertCount
	header maxLinkCount = maxLinkCount
	vcopy(header bmin, params bmin)
	vcopy(header bmax, params bmax)
	header detailMeshCount = params polyCount
	header detailVertCount = uniqueDetailVertCount
	header detailTriCount = params detailTriCount
	header bvQuantFactor = 1.0 / params cs
	header offMeshBase = params polyCount
	header walkableHeight = params walkableHeight
	header walkableRadius = params walkableRadius
	header walkableClimb = params walkableClimb
	header offMeshConCount = storedOffMeshConCount
	header bvNodeCount = params polyCount*2
	
	offMeshVertsBase = params vertCount
	offMeshPolyBase = params polyCount
	
	// Store vertices
	// Mesh vertices
	for (Int i = 0; i < params vertCount; ++i) {
		iv: UShort* = &params verts[i*3]
		v: Float* = &navVerts[i*3]
		v[0] = params bmin[0] + iv[0] * params cs
		v[1] = params bmin[1] + iv[1] * params ch
		v[2] = params bmin[2] + iv[2] * params cs
	}
	// Off-mesh link vertices.
	n := 0
	for (i: Int in 0..(params offMeshConCount)) {
		// Only store connections which start from this tile.
		if (offMeshConFlags[i*2+0] == 0xff) {
			linkv: Float* = &params offMeshConVerts[i*2*3]
			v: Float* = &navVerts[(offMeshVertsBase + n*2)*3]
			vcopy(&v[0], &linkv[0])
			vcopy(&v[3], &linkv[3])
			n++
		}
	}
	
	// Store polygons
	// Mesh polys
	src: UShort* = params polys
	for (i: Int in 0..(params polyCount)) {
		p: DTPoly* = &navPolys[i]
		p vertCount = 0
		p flags = DT_POLY_GROUND
		for (j: Int in 0..nvp) {
			if (src[j] == MESH_NULL_IDX) break
			p verts[j] = src[j]
			p neis[j] = (src[nvp+j]+1) & 0xffff
			p vertCount++
		}
		src += nvp*2
	}
	// Off-mesh connection vertices.
	n = 0
	for (i: Int in 0..(params offMeshConCount)) {
		// Only store connections which start from this tile.
		if (offMeshConFlags[i*2+0] == 0xff) {
			p: DTPoly* = &navPolys[offMeshPolyBase+n]
			p vertCount = 2
			p verts[0] = (offMeshVertsBase + n*2+0) as UShort
			p verts[1] = (offMeshVertsBase + n*2+1) as UShort
			p flags = DT_POLY_OFFMESH_CONNECTION; // Off-mesh link poly.
			n++
		}
	}
	
	// Store portal edges.
	if (params tileSize > 0) {
		for (i: Int in 0..(params polyCount)) {
			poly: DTPoly* = &navPolys[i]
			for (j: Int in 0..(poly vertCount)) {
				nj: Int = j+1
				if (nj >= poly vertCount) nj = 0
				va: UShort* = &params verts[poly verts[j]*3]
				vb: UShort* = &params verts[poly verts[nj]*3]
				
				if (va[0] == params tileSize && vb[0] == params tileSize) // x+
					poly neis[j] = DT_EXT_LINK | 0
				else if (va[2] == params tileSize && vb[2]  == params tileSize) // z+
					poly neis[j] = DT_EXT_LINK | 2
				else if (va[0] == 0 && vb[0] == 0) // x-
					poly neis[j] = DT_EXT_LINK | 4
				else if (va[2] == 0 && vb[2] == 0) // z-
					poly neis[j] = DT_EXT_LINK | 6
			}
		}
	}

	// Store detail meshes and vertices.
	// The nav polygon vertices are stored as the first vertices on each mesh.
	// We compress the mesh data by skipping them and using the navmesh coordinates.
	vbase: UShort = 0
	for (i: Int in 0..(params polyCount)) {
		dtl: DTPolyDetail& = navDMeshes[i]
		nb: Int = params detailMeshes[i*4+0]
		ndv := params detailMeshes[i*4+1]
		nv: Int = navPolys[i] vertCount
		dtl vertBase = vbase
		dtl vertCount = ndv-nv
		dtl triBase = params detailMeshes[i*4+2]
		dtl triCount = params detailMeshes[i*4+3]
		// Copy vertices except the first 'nv' verts which are equal to nav poly verts.
		if (ndv-nv) {
			memcpy(&navDVerts[vbase*3], &params detailVerts[(vb+nv)*3], sizeof(Float)*3*(ndv-nv))
			vbase += ndv-nv
		}
	}
	// Store triangles.
	memcpy(navDTris, params detailTris, sizeof(UInt8)*4*params detailTriCount)
	
	// Store and create BVtree.
	// TODO: take detail mesh into account! use byte per bbox extent?
	createBVTree(params verts, params vertCount, params polys, params polyCount,
				 nvp, params cs, params ch, params polyCount*2, navBvtree)
	
	// Store Off-Mesh connections.
	n = 0
	for (i: Int in 0..(params offMeshConCount)) {
		// Only store connections which start from this tile.
		if (offMeshConFlags[i*2+0] == 0xff) {
			con: DTOffMeshConnection* = &offMeshCons[n]
			con poly = offMeshPolyBase + n
			// Copy connection end-points.
			endPts: Float* = &params offMeshConVerts[i*2*3]
			vcopy(&con pos[0], &endPts[0])
			vcopy(&con pos[3], &endPts[3])
			con rad = params offMeshConRad[i]
			con flags = params offMeshConDir[i] ? DT_OFFMESH_CON_BIDIR : 0
			con side = offMeshConFlags[i*2+1]
			n++
		}
	}
	
	delete [] offMeshConFlags
	
	*outData = data
	*outDataSize = dataSize
	return true
}

swapByte: inline func(a, b: UInt8*) {
	tmp: UInt8 = *a
	*a = *b
	*b = tmp
}

swapEndian: inline func(v: UShort*) {
	x := v as UInt8*
	swapByte(x+0, x+1)
}

swapEndian: inline func(v: Short*) {
	x := v as UInt8*
	swapByte(x+0, x+1)
}

swapEndian: inline func(UInt* v) {
	x := v as UInt8*
	swapByte(x+0, x+3); swapByte(x+1, x+2)
}

swapEndian: inline func(Int* v) {
	x := v as UInt8*
	swapByte(x+0, x+3); swapByte(x+1, x+2)
}

swapEndian: inline func(Float* v) {
	x := v as UInt8*
	swapByte(x+0, x+3); swapByte(x+1, x+2)
}

dtNavMeshHeaderSwapEndian: func(data: UInt8*, dataSize: Int) -> Bool {
	header: DTMeshHeader* = data as DTMeshHeader*
	
	swappedMagic := DT_NAVMESH_MAGIC
	swappedVersion := DT_NAVMESH_VERSION
	swapEndian(&swappedMagic)
	swapEndian(&swappedVersion)
	
	if ((header magic != DT_NAVMESH_MAGIC || header dversion != DT_NAVMESH_VERSION) &&
		(header magic != swappedMagic || header dversion != swappedVersion)) {
		return false
	}
	
	swapEndian(&header magic)
	swapEndian(&header dversion)
	swapEndian(&header polyCount)
	swapEndian(&header vertCount)
	swapEndian(&header maxLinkCount)
	swapEndian(&header detailMeshCount)
	swapEndian(&header detailVertCount)
	swapEndian(&header detailTriCount)
	swapEndian(&header bvNodeCount)
	swapEndian(&header offMeshConCount)
	swapEndian(&header offMeshBase)
	swapEndian(&header walkableHeight)
	swapEndian(&header walkableRadius)
	swapEndian(&header walkableClimb)
	swapEndian(&header bmin[0])
	swapEndian(&header bmin[1])
	swapEndian(&header bmin[2])
	swapEndian(&header bmax[0])
	swapEndian(&header bmax[1])
	swapEndian(&header bmax[2])
	swapEndian(&header bvQuantFactor)
	
	// Freelist index and pointers are updated when tile is added, no need to swap.
	return true
}

dtNavMeshDataSwapEndian: func(UInt8* data, Int dataSize) -> Bool {
	// Make sure the data is in right format.
	header: DTMeshHeader* = data as DTMeshHeader*
	if (header magic != DT_NAVMESH_MAGIC)
		return false
	if (header dversion != DT_NAVMESH_VERSION)
		return false
	
	// Patch header pointers.
	headerSize := align4(sizeof(DTMeshHeader))
	vertsSize := align4(sizeof(Float)*3*header vertCount)
	polysSize := align4(sizeof(DTPoly)*header polyCount)
	linksSize := align4(sizeof(DTLink)*(header maxLinkCount))
	detailMeshesSize := align4(sizeof(DTPolyDetail)*header detailMeshCount)
	detailVertsSize := align4(sizeof(Float)*3*header detailVertCount)
	detailTrisSize := align4(sizeof(UInt8)*4*header detailTriCount)
	bvtreeSize := align4(sizeof(DTBVNode)*header bvNodeCount)
	offMeshLinksSize := align4(sizeof(DTOffMeshConnection)*header offMeshConCount)
	
	d: UInt8* = data + headerSize
	verts: Float* = d as Float*; d += vertsSize
	polys: DTPoly* = d as DTPoly*; d += polysSize
	/*links: DTLink* = d as DTLink*);*/ d += linksSize
	detailMeshes: DTPolyDetail* = d as DTPolyDetail*; d += detailMeshesSize
	detailVerts: Float* = d as Float*; d += detailVertsSize
	/*detailTris: UInt8* = d as UInt8*;*/ d += detailTrisSize
	bvTree: DTBVNode* = d as DTBVNode*; d += bvtreeSize
	offMeshCons: DTOffMeshConnection* = d as DTOffMeshConnection*; d += offMeshLinksSize
	
	// Vertices
	for (i: Int in 0..(header vertCount*3)) {
		swapEndian(&verts[i])
	}
	
	// Polys
	for (i: Int in 0..(header polyCount)) {
		p: DTPoly* = &polys[i]
		// poly firstLink is update when tile is added, no need to swap.
		for (Int j = 0; j < DT_VERTS_PER_POLYGON; ++j) {
			swapEndian(&p verts[j])
			swapEndian(&p neis[j])
		}
		swapEndian(&p flags)
	}
	
	// Links are rebuild when tile is added, no need to swap.
	// Detail meshes
	for (i: Int in 0..(header detailMeshCount)) {
		DTPolyDetail* pd = &detailMeshes[i]
		swapEndian(&pd vertBase)
		swapEndian(&pd vertCount)
		swapEndian(&pd triBase)
		swapEndian(&pd triCount)
	}
	
	// Detail verts
	for (i: Int in 0..(header detailVertCount*3)) {
		swapEndian(&detailVerts[i])
	}
	
	// BV-tree
	for (i: Int in 0..(header bvNodeCount)) {
		DTBVNode* node = &bvTree[i]
		for (Int j = 0; j < 3; ++j) {
			swapEndian(&node bmin[j])
			swapEndian(&node bmax[j])
		}
		swapEndian(&node i)
	}
	
	// Off-mesh Connections.
	for (i: Int in 0..(header offMeshConCount)) {
		con: DTOffMeshConnection* = &offMeshCons[i]
		for (Int j = 0; j < 6; ++j)
			swapEndian(&con pos[j])
		swapEndian(&con rad)
		swapEndian(&con poly)
	}
	
	return true
}

