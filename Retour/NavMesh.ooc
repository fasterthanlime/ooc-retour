
opposite: inline func(side: Int) -> Int {
	return (side+4) & 0x7
}

overlapBoxes: inline func(amin, amax, bmin, bmax: Float*) -> Bool {
	overlap := true
	overlap = (amin[0] > bmax[0] || amax[0] < bmin[0]) ? false : overlap
	overlap = (amin[1] > bmax[1] || amax[1] < bmin[1]) ? false : overlap
	overlap = (amin[2] > bmax[2] || amax[2] < bmin[2]) ? false : overlap
	return overlap
}

overlapRects: inline func(amin, amax, bmin, bmax: Float*) -> Bool {
	overlap := true
	overlap = (amin[0] > bmax[0] || amax[0] < bmin[0]) ? false : overlap
	overlap = (amin[1] > bmax[1] || amax[1] < bmin[1]) ? false : overlap
	return overlap
}

calcRect: static func(va, vb, bmin, bmax: Float*, side: Int, padx, pady: Float) {
	if (side == 0 || side == 4) {
		bmin[0] = min(va[2], vb[2]) + padx
		bmin[1] = min(va[1], vb[1]) - pady
		bmax[0] = max(va[2], vb[2]) - padx
		bmax[1] = max(va[1], vb[1]) + pady
	} else if (side == 2 || side == 6) {
		bmin[0] = min(va[0], vb[0]) + padx
		bmin[1] = min(va[1], vb[1]) - pady
		bmax[0] = max(va[0], vb[0]) - padx
		bmax[1] = max(va[1], vb[1]) + pady
	}
}

computeTileHash: inline func(x, y, mask: Int) {
	h1: static const UInt = 0x8da6b343 // Large multiplicative constants
	h2: static const UInt = 0xd8163841 // here arbitrarily chosen primes
	n: UInt = h1 * x + h2 * y
	return (n & mask) as Int
}

allocLink inline func(tile: DTMeshTile*) -> UInt {
	if (tile header linksFreeList == DT_NULL_LINK)
		return DT_NULL_LINK
	link := tile header linksFreeList
	tile header linksFreeList = tile header links[link] next
	return link
}

freeLink: inline func(tile: DTMeshTile*, link: UInt) {
	tile header links[link] next = tile header linksFreeList
	tile header linksFreeList = link
}

passFilter: inline func(filter: DTQueryFilter*, flags: UShort) -> Bool {
	return (flags & filter includeFlags) != 0 && (flags & filter excludeFlags) == 0
}

// Reference to navigation polygon.
//typedef UInt DTPolyRef
DTPolyRef: cover from UInt // correcto?

// Maximum number of vertices per navigation polygon.
DT_VERTS_PER_POLYGON: static const Int = 6

DT_NAVMESH_MAGIC: static const Int = 'DNAV' // huh?
DT_NAVMESH_VERSION: static const Int = 2

DT_EXT_LINK: static const UShort = 0x8000
DT_NULL_LINK: static const UInt = 0xffffffff
DT_OFFMESH_CON_BIDIR: static const UInt = 1

// Flags returned by findStraightPath().
//enum
DTStraightPathFlags: class {
	DT_STRAIGHTPATH_START: static const Int = 0x01				// The vertex is the start position.
	DT_STRAIGHTPATH_END: static const Int = 0x02					// The vertex is the end position.
	DT_STRAIGHTPATH_OFFMESH_CONNECTION: static const Int = 0x04	// The vertex is start of an off-mesh link.
}

// Flags describing polygon properties.
//enum
DTPolyFlags: class {
	DT_POLY_GROUND: static const Int = 0x01						// Regular ground polygons.
	DT_POLY_OFFMESH_CONNECTION: static const Int = 0x02			// Off-mesh connections.
}

DTQueryFilter: class {
	init: func() {
		includeFlags = 0xffff
		excludeFlags = 0
	}
	
	includeFlags: UShort		// If any of the flags are set, the poly is included.
	excludeFlags: UShort		// If any of the flags are set, the poly is excluded.
}

// Structure describing the navigation polygon data.
DTPoly: class {
	firstLink: UInt							// Index to first link in linked list. 
	verts: UShort[DT_VERTS_PER_POLYGON]		// Indices to vertices of the poly.
	neis: UShort[DT_VERTS_PER_POLYGON]		// Refs to neighbours of the poly.
	flags: UShort					// Flags (see DTPolyFlags).
	vertCount: UInt8				// Number of vertices.
}

// Stucture describing polygon detail triangles.
DTPolyDetail: class {
	vertBase: UShort		// Offset to detail vertex array.
	vertCount: UShort		// Number of vertices in the detail mesh.
	triBase: UShort			// Offset to detail triangle array.
	triCount: UShort		// Number of triangles.
}

// Stucture describing a link to another polygon.
DTLink: class {
	ref: DTPolyRef				// Neighbour reference.
	next: UInt					// Index to next link.
	edge: UInt8					// Index to polygon edge which owns this link. 
	side: UInt8					// If boundary link, defines on which side the link is.
	bmin, bmax: UInt8			// If boundary link, defines the sub edge area.
}

DTBVNode: class {
	bmin, bmax: UInt8[3]		// BVnode bounds
	i: Int						// Index to item or if negative, escape index.
}

DTOffMeshConnection: class {
	pos: Float[6]			// Both end point locations.
	rad: Float				// Link connection radius.
	poly: UShort			// Poly Id
	flags: UInt8			// Link flags
	side: UInt8				// End point side.
}

DTMeshHeader: class {
	magic: Int								// Magic number, used to identify the data.
	dversion: Int							// Data version number.
	
	polyCount: Int							// Number of polygons in the tile.
	vertCount: Int							// Number of vertices in the tile.
	maxLinkCount: Int						// Number of allocated links.
	detailMeshCount: Int					// Number of detail meshes.
	detailVertCount: Int					// Number of detail vertices.
	detailTriCount: Int						// Number of detail triangles.
	bvNodeCount: Int						// Number of BVtree nodes.
	offMeshConCount: Int					// Number of Off-Mesh links.
	offMeshBase: Int						// Index to first polygon which is Off-Mesh link.
	linksFreeList: UInt						// Index to next free link.
	walkableHeight: Float					// Height of the agent.
	walkableRadius: Float					// Radius of the agent
	walkableClimb: Float					// Max climb height of the agent.
	bmin, bmax: Float[3]					// Bounding box of the tile.
	bvQuantFactor: Float					// BVtree quantization factor (world to bvnode coords)
	polys: DTPoly*							// Pointer to the polygons (will be updated when tile is added).
	verts: Float*							// Pointer to the vertices (will be updated when tile added).
	links: DTLink*							// Pointer to the links (will be updated when tile added).
	detailMeshes: DTPolyDetail*				// Pointer to detail meshes (will be updated when tile added).
	detailVerts: Float*						// Pointer to detail vertices (will be updated when tile added).
	detailTris: UInt8*						// Pointer to detail triangles (will be updated when tile added).
	bvTree: DTBVNode*						// Pointer to BVtree nodes (will be updated when tile added).
	offMeshCons: DTOffMeshConnection*		// Pointer to Off-Mesh links. (will be updated when tile added).
}

DTMeshTile: class {
	salt: Int								// Counter describing modifications to the tile.
	x, y: Int								// Grid location of the tile.
	header: DTMeshHeader*					// Pointer to tile header.
	data: UInt8*					// Pointer to tile data.
	dataSize: Int							// Size of the tile data.
	ownsData: Bool							// Flag indicating of the navmesh should release the data.
	next: DTMeshTile*						// Next free tile or, next tile in spatial grid.
}

DTNavMesh: class {
	
	init: func() {
		m_orig[0] = 0
		m_orig[1] = 0
		m_orig[2] = 0
	}
	
	destroy: func() {
		for (i: Int in 0..m_maxTiles) {
			if (m_tiles[i] data && m_tiles[i] ownsData) {
				delete [] m_tiles[i] data
				m_tiles[i] data = 0
				m_tiles[i] dataSize = 0
			}
		}
		delete m_nodePool
		delete m_openList
		delete [] m_posLookup
		delete [] m_tiles
	}
	
	initialize: func~multiTile(orig: Float*, tileWidth, tileHeight: Float, maxTiles, maxPolys, maxNodes: Int) -> Bool {
		vcopy(m_orig, orig)
		m_tileWidth = tileWidth
		m_tileHeight = tileHeight
		
		// Init tiles
		m_maxTiles = maxTiles
		m_tileLutSize = nextPow2(maxTiles/4)
		if (!m_tileLutSize)
			m_tileLutSize = 1
		m_tileLutMask = m_tileLutSize-1
		
		m_tiles = rcAllocArray(DTMeshTile, m_maxTiles)
		if (!m_tiles)
			return false
		
		m_posLookup = rcAllocArray(DTMeshTile*, m_tileLutSize)
		if (!m_posLookup)
			return false
		
		memset(m_tiles, 0, sizeof(DTMeshTile)*m_maxTiles)
		memset(m_posLookup, 0, sizeof(DTMeshTile*)*m_tileLutSize)
		m_nextFree = 0
		for (i := m_maxTiles-1; i >= 0; --i) {
			m_tiles[i] next = m_nextFree
			m_nextFree = &m_tiles[i]
		}
		
		if (!m_nodePool) {
			m_nodePool = new DTNodePool(maxNodes, nextPow2(maxNodes/4))
			if (!m_nodePool)
				return false
		}
		
		if (!m_openList) {
			m_openList = new DTNodeQueue(maxNodes)
			if (!m_openList)
				return false
		}
		
		// Init ID generator values.
		m_tileBits = max(1 as UInt, ilog2(nextPow2(maxTiles as UInt)))
		m_polyBits = max(1 as UInt, ilog2(nextPow2(maxPolys as UInt)))
		m_saltBits = 32 - m_tileBits - m_polyBits
		if (m_saltBits < 10)
			return false
		
		return true
	}
	
	initialize: func~singleTile(data: UInt8*, dataSize: Int, ownsData: Bool, maxNodes: Int) -> Bool {
		// Make sure the data is in right format.
		header: DTMeshHeader* = (DTMeshHeader*)data
		if (header magic != DT_NAVMESH_MAGIC)
			return false
		if (header dversion != DT_NAVMESH_VERSION)
			return false
		
		w: Float = header bmax[0] - header bmin[0]
		h: Float = header bmax[2] - header bmin[2]
		if (!init(header bmin, w, h, 1, header polyCount, maxNodes))
			return false

		return addTileAt(0, 0, data, dataSize, ownsData)
	}
	
	findConnectingPolys: func(va, vb: Float*, tile: DTMeshTile*, side: Int, con: DTPolyRef*, conarea: Float*, maxcon: Int) -> Int {
		if (!tile) return 0
		h: DTMeshHeader* = tile header
		amin, amax: Float[2]
		calcRect(va, vb, amin, amax, side, 0.01, h walkableClimb)
		
		// Remove links pointing to 'side' and compact the links array. 
		bmin, bmax: Float[2]
		m: UShort = DT_EXT_LINK | (side as UShort)
		n := 0
		
		base := getTileId(tile)
		for (i: Int in 0..(h polyCount)) {
			poly: DTPoly* = &h polys[i]
			nv := poly vertCount
			for (j: Int in 0..nv) {
				// Skip edges which do not point to the right side.
				if (poly neis[j] != m) continue
				// Check if the segments touch.
				vc: Float* = &h verts[poly verts[j]*3]
				vd: Float* = &h verts[poly verts[(j+1) % nv]*3]
				calcRect(vc, vd, bmin, bmax, side, 0.01, h walkableClimb)
				if (!overlapRects(amin, amax, bmin, bmax)) continue
				// Add return value.
				if (n < maxcon) {
					conarea[n*2+0] = max(amin[0], bmin[0])
					conarea[n*2+1] = min(amax[0], bmax[0])
					con[n] = base | (i as UInt)
					n++
				}
				break
			}
		}
		return n
	}
	
	unconnectExtLinks: func(tile: DTMeshTile*, side: Int) {
		if (!tile)
			return
		
		header: DTMeshHeader* = tile header
		for (i: Int in 0..(header polyCount)) {
			poly: DTPoly* = &header polys[i]
			j := poly firstLink
			pj := DT_NULL_LINK
			while (j != DT_NULL_LINK) {
				if (header links[j] side == side) {
					// Revove link.
					nj := header links[j] next
					if (pj == DT_NULL_LINK)
						poly firstLink = nj
					else
						header links[pj] next = nj
					freeLink(tile, j)
					j = nj
				} else {
					// Advance
					pj = j
					j = header links[j] next
				}
			}
		}
	}
	
	connectExtLinks: func(tile: DTMeshTile*, target: DTMeshTile*, side: Int) {
		if (!tile)
			return
		
		header: DTMeshHeader* = tile header
		// Connect border links.
		for (i: Int in 0..(header polyCount)) {
			poly: DTPoly* = &header polys[i]
			
			// Create new links.
			m: UShort = DT_EXT_LINK | (side as UShort)
			nv := poly vertCount
			for (j: Int in 0..nv) {
				// Skip edges which do not point to the right side.
				if (poly neis[j] != m)
					continue
				
				// Create new links
				va: Float* = &header verts[poly verts[j]*3]
				vb: Float* = &header verts[poly verts[(j+1) % nv]*3]
				nei: DTPolyRef[4]
				neia: Float[4*2]
				nnei := findConnectingPolys(va, vb, target, opposite(side), nei, neia, 4)
				for (k: Int in 0..nnei) {
					idx := allocLink(tile)
					if (idx != DT_NULL_LINK) {
						link: DTLink* = &header links[idx]
						link ref = nei[k]
						link edge = j as UInt8
						link side = side as UInt8
						
						link next = poly firstLink
						poly firstLink = idx
						
						// Compress portal limits to a byte value.
						if (side == 0 || side == 4) {
							lmin := min(va[2], vb[2])
							lmax := max(va[2], vb[2])
							link bmin = (clamp((neia[k*2+0]-lmin)/(lmax-lmin), 0.0, 1.0)*255.0) as UInt8
							link bmax = (clamp((neia[k*2+1]-lmin)/(lmax-lmin), 0.0, 1.0)*255.0) as UInt8
						} else if (side == 2 || side == 6) {
							lmin := min(va[0], vb[0])
							lmax := max(va[0], vb[0])
							link bmin = (clamp((neia[k*2+0]-lmin)/(lmax-lmin), 0.0, 1.0)*255.0) as UInt8
							link bmax = (clamp((neia[k*2+1]-lmin)/(lmax-lmin), 0.0, 1.0)*255.0) as UInt8
						}
					}
				}
			}
		}
	}

	connectExtOffMeshLinks: func(tile, target: DTMeshTile*, side: Int) {
		if (!tile)
			return
		
		header: DTMeshHeader* = tile header
		// Connect off-mesh links.
		// We are interested on links which land from target tile to this tile.
		targetHeader: DTMeshHeader* = target header
		oppositeSide := opposite(side) as UInt8
		defaultFilter: DTQueryFilter
		
		for (i: Int in 0..(targetHeader offMeshConCount)) {
			DTOffMeshConnection* targetCon = &targetHeader offMeshCons[i]
			if (targetCon side != oppositeSide)
				continue
			
			targetPoly: DTPoly* = &targetHeader polys[targetCon poly]
			ext: Float[3] = [targetCon rad, targetHeader walkableClimb, targetCon rad]
			
			// Find polygon to connect to.
			p: Float* = &targetCon pos[3]
			nearestPt: Float[3]
			ref := findNearestPolyInTile(tile, p, ext, &defaultFilter, nearestPt)
			if (!ref) continue
			// findNearestPoly may return too optimistic results, further check to make sure. 
			if (sqr(nearestPt[0]-p[0])+sqr(nearestPt[2]-p[2]) > sqr(targetCon rad))
				continue
			// Make sure the location is on current mesh.
			v: Float* = &targetHeader verts[targetPoly verts[1]*3]
			vcopy(v, nearestPt)
			
			// Link off-mesh connection to target poly.
			idx := allocLink(target)
			if (idx != DT_NULL_LINK) {
				link: DTLink* = &targetHeader links[idx]
				link ref = ref
				link edge = 1 as UInt8
				link side = oppositeSide
				link bmin = link bmax = 0
				// Add to linked list.
				link next = targetPoly firstLink
				targetPoly firstLink = idx
			}
		
			// Link target poly to off-mesh connection.
			if (targetCon flags & DT_OFFMESH_CON_BIDIR) {
				idx: UInt = allocLink(tile)
				if (idx != DT_NULL_LINK) {
					landPolyIdx: UShort = decodePolyIdPoly(ref)
					landPoly: DTPoly* = &header polys[landPolyIdx]
					link: DTLink* = &header links[idx]
					link ref = getTileId(target) | (targetCon poly as UInt)
					link edge = 0
					link side = side
					link bmin = link bmax = 0
					// Add to linked list.
					link next = landPoly firstLink
					landPoly firstLink = idx
				}
			}
		}

	}
	
	connectIntLinks: func(tile: DTMeshTile*) {
		if (!tile)
			return
		
		header: DTMeshHeader* = tile header
		base: DTPolyRef = getTileId(tile)
		for (i: Int in 0..(header polyCount)) {
			poly: DTPoly* = &header polys[i]
			poly firstLink = DT_NULL_LINK
			
			if (poly flags & DT_POLY_OFFMESH_CONNECTION)
				continue
			
			// Build edge links backwards so that the links will be
			// in the linked list from lowest index to highest.
			for (j: Int = (poly vertCount-1); j >= 0; --j) {
				// Skip hard and non-internal edges.
				if (poly neis[j] == 0 || (poly neis[j] & DT_EXT_LINK)) continue
				
				idx: UInt = allocLink(tile)
				if (idx != DT_NULL_LINK) {
					link: DTLink* = &header links[idx]
					link ref = base | (UInt)(poly neis[j]-1)
					link edge = (UInt8)j
					link side = 0xff
					link bmin = link bmax = 0
					// Add to linked list.
					link next = poly firstLink
					poly firstLink = idx
				}
			}			
		}
	}
	
	connectIntOffMeshLinks: func(tile: DTMeshTile*) {
		if (!tile) return
		header: DTMeshHeader* = tile header
		base: DTPolyRef = getTileId(tile)
		
		// Find Off-mesh connection end points.
		for (i: Int in 0..(header offMeshConCount)) {
			con: DTOffMeshConnection* = &header offMeshCons[i]
			poly: DTPoly* = &header polys[con poly]
			defaultFilter: DTQueryFilter
			
			ext: Float[3] = { con rad, header walkableClimb, con rad }
			
			for (j: Int in 0..2) {
				side: UInt8 = j == 0 ? 0xff : con side
				
				if (side == 0xff) {
					// Find polygon to connect to.
					p: Float* = &con pos[j*3]
					nearestPt: Float[3]
					ref: DTPolyRef = findNearestPolyInTile(tile, p, ext, &defaultFilter, nearestPt)
					if (!ref) continue
					// findNearestPoly may return too optimistic results, further check to make sure. 
					if (sqr(nearestPt[0]-p[0])+sqr(nearestPt[2]-p[2]) > sqr(con rad))
						continue
					// Make sure the location is on current mesh.
					v: Float* = &header verts[poly verts[j]*3]
					vcopy(v, nearestPt)
					
					// Link off-mesh connection to target poly.
					idx: UInt = allocLink(tile)
					if (idx != DT_NULL_LINK) {
						link: DTLink* = &header links[idx]
						link ref = ref
						link edge = j as UInt8
						link side = 0xff
						link bmin = link bmax = 0
						// Add to linked list.
						link next = poly firstLink
						poly firstLink = idx
					}
					
					// Start end-point is always connect back to off-mesh connection, 
					// Destination end-point only if it is bidirectional link. 
					if (j == 0 || (j == 1 && (con flags & DT_OFFMESH_CON_BIDIR))) {
						// Link target poly to off-mesh connection.
						idx: UInt = allocLink(tile)
						if (idx != DT_NULL_LINK) {
							landPolyIdx: UShort = decodePolyIdPoly(ref)
							landPoly: DTPoly* = &header polys[landPolyIdx]
							link: DTLink* = &header links[idx]
							link ref = base | (UInt)(con poly)
							link edge = 0
							link side = 0xff
							link bmin = link bmax = 0
							// Add to linked list.
							link next = landPoly firstLink
							landPoly firstLink = idx
						}
					}
				
				}
			}
		}
	}
	
	addTileAt: func(x, y: Int, data: UInt8*, dataSize: Int, ownsData: Bool) -> Bool {
		if (getTileAt(x, y))
			return false
		// Make sure there is enough space for new tile.
		if (!m_nextFree)
			return false
		// Make sure the data is in right format.
		header: DTMeshHeader* = (DTMeshHeader*)data
		if (header magic != DT_NAVMESH_MAGIC)
			return false
		if (header dversion != DT_NAVMESH_VERSION)
			return false
		
		// Allocate a tile.
		tile: DTMeshTile* = m_nextFree
		m_nextFree = tile next
		tile next = 0
		
		// Insert tile into the position lut.
		h := computeTileHash(x, y, m_tileLutMask)
		tile next = m_posLookup[h]
		m_posLookup[h] = tile
	
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
		header verts = d as Float*; d += vertsSize
		header polys = d as DTPoly*; d += polysSize
		header links = d as DTLink*; d += linksSize
		header detailMeshes = d as DTPolyDetail*; d += detailMeshesSize
		header detailVerts = d as Float*; d += detailVertsSize
		header detailTris = d as UInt8*; d += detailTrisSize
		header bvTree = d as DTBVNode*; d += bvtreeSize
		header offMeshCons = d as DTOffMeshConnection*; d += offMeshLinksSize
		
		// Build links freelist
		header linksFreeList = 0
		header links[header maxLinkCount-1] next = DT_NULL_LINK
		for (i: Int in 0..(header maxLinkCount-1))
			header links[i] next = i+1
		
		// Init tile.
		tile header = header
		tile x = x
		tile y = y
		tile data = data
		tile dataSize = dataSize
		tile ownsData = ownsData

		connectIntLinks(tile)
		connectIntOffMeshLinks(tile)

		// Create connections connections.
		for (i: Int in 0..8) {
			nei: DTMeshTile* = getNeighbourTileAt(x, y, i)
			if (nei) {
				connectExtLinks(tile, nei, i)
				connectExtLinks(nei, tile, opposite(i))
				connectExtOffMeshLinks(tile, nei, i)
				connectExtOffMeshLinks(nei, tile, opposite(i))
			}
		}
		return true
	}

	getTileAt: func(x, y: Int) -> DTMeshTile* {
		// Find tile based on hash.
		h := computeTileHash(x, y, m_tileLutMask)
		tile: DTMeshTile* = m_posLookup[h]
		while (tile) {
			if (tile x == x && tile y == y)
				return tile
			tile = tile next
		}
		return 0
	}

	getMaxTiles: func() -> Int {
		return m_maxTiles
	}

	getTile: func(Int i) -> DTMeshTile* {
		return &m_tiles[i]
	}

	getTile: func(Int i) -> DTMeshTile* {
		return &m_tiles[i]
	}
	
	getTileByRef: func(ref: DTPolyRef, polyIndex: Int*) -> DTMeshTile* {
		salt, it, ip: UInt
		decodePolyId(ref, salt, it, ip)
		if (it >= m_maxTiles as UInt) return 0
		if (m_tiles[it] salt != salt || m_tiles[it] header == 0) return 0
		if (ip >= (m_tiles[it] header polyCount as UInt)) return 0
		if (polyIndex) *polyIndex = (Int)ip
		return &m_tiles[it]
	}

	getNeighbourTileAt: func(x, y, side: Int) -> DTMeshTile* {
		match (side) {
			case 0 => x++; break
			case 1 => x++; y++; break
			case 2 => y++; break
			case 3 => x--; y++; break
			case 4 => x--; break
			case 5 => x--; y--; break
			case 6 => y--; break
			case 7 => x++; y--; break
		}
		return getTileAt(x, y)
	}

	removeTileAt: func(x, y: Int, data: UInt8**, dataSize: Int*) -> Bool {
		// Remove tile from hash lookup.
		h := computeTileHash(x, y, m_tileLutMask)
		prev: DTMeshTile* = 0
		tile: DTMeshTile* = m_posLookup[h]
		while (tile) {
			if (tile x == x && tile y == y) {
				if (prev)
					prev next = tile next
				else
					m_posLookup[h] = tile next
				break
			}
			prev = tile
			tile = tile next
		}
		if (!tile)
			return false
		
		// Remove connections to neighbour tiles.
		for (i: Int in 0..8) {
			nei: DTMeshTile* = getNeighbourTileAt(x, y, i)
			if (!nei) continue
			unconnectExtLinks(nei, opposite(i))
		}
		
		// Reset tile.
		if (tile ownsData) {
			// Owns data
			delete [] tile data
			tile data = 0
			tile dataSize = 0
			if (data) *data = 0
			if (dataSize) *dataSize = 0
		} else {
			if (data) *data = tile data
			if (dataSize) *dataSize = tile dataSize
		}
		tile header = 0
		tile x = tile y = 0
		tile salt++

		// Add to free list.
		tile next = m_nextFree
		m_nextFree = tile

		return true
	}

	getTileId: func(tile: DTMeshTile*) -> DTPolyRef{
		if (!tile) return 0
		it: UInt = tile - m_tiles
		return encodePolyId(tile salt, it, 0)
	}
	
	closestPointOnPoly: func(ref: DTPolyRef, pos, closest: Float*) -> Bool {
		salt, it, ip: UInt
		decodePolyId(ref, salt, it, ip)
		if (it >= m_maxTiles as UInt) return false
		if (m_tiles[it] salt != salt || m_tiles[it] header == 0) return false
		header: DTMeshHeader* = m_tiles[it] header
		if (ip >= (header polyCount as UInt)) return false
		
		return closestPointOnPolyInTile(&m_tiles[it], ip, pos, closest)
	}

	closestPointOnPolyInTile: func(tile: DTMeshTile*, ip: UInt, pos, closest: Float*) -> Bool {
		header: DTMeshHeader* = tile header
		poly: DTPoly* = &header polys[ip]
		
		closestDistSqr := FLT_MAX
		pd: DTPolyDetail* = &header detailMeshes[ip]
		
		for (j: Int in 0..(pd triCount)) {
			t: UInt8* = &header detailTris[(pd triBase+j)*4]
			v: Float*[3]
			for (k: Int in 0..3) {
				if (t[k] < poly vertCount)
					v[k] = &header verts[poly verts[t[k]]*3]
				else
					v[k] = &header detailVerts[(pd vertBase+(t[k]-poly vertCount))*3]
			}
			p: Float[3]
			closestPtPointTriangle(pt, pos, v[0], v[1], v[2])
			d: Float = vdistSqr(pos, pt)
			if (d < closestDistSqr) {
				vcopy(closest, pt)
				closestDistSqr = d
			}
		}
		return true
	}
	
	closestPointOnPolyBoundary: func(ref: DTPolyRef, pos, closest: Float*) -> Bool {
		salt, it, ip: UInt
		decodePolyId(ref, salt, it, ip)
		if (it >= m_maxTiles as UInt) return false
		if (m_tiles[it] salt != salt || m_tiles[it] header == 0) return false
		header: DTMeshHeader* = m_tiles[it] header
		
		if (ip >= (header polyCount as UInt)) return false
		poly: DTPoly* = &header polys[ip]
		
		// Collect vertices.
		verts: Float[DT_VERTS_PER_POLYGON*3]
		edged: Float[DT_VERTS_PER_POLYGON]
		edget: Float[DT_VERTS_PER_POLYGON]
		nv := 0
		for (i: Int in 0..(poly vertCount as Int)) {
			vcopy(&verts[nv*3], &header verts[poly verts[i]*3])
			nv++
		}		
		
		inside := distancePtPolyEdgesSqr(pos, verts, nv, edged, edget)
		if (inside) {
			// Point is inside the polygon, return the point.
			vcopy(closest, pos)
		} else {
			// Point is outside the polygon, clamp to nearest edge.
			dmin := FLT_MAX
			imin := -1
			for (i: Int in 0..nv) {
				if (edged[i] < dmin) {
					dmin = edged[i]
					imin = i
				}
			}
			va: Float* = &verts[imin*3]
			vb: Float* = &verts[((imin+1)%nv)*3]
			vlerp(closest, va, vb, edget[imin])
		}
		return true
	}

	// Returns start and end location of an off-mesh link polygon.
	getOffMeshConnectionPolyEndPoints: func(prevRef, polyRef: DTPolyRef, startPos, endPos: Float*) -> Bool {
		salt, it, ip: UInt
		
		// Get current polygon
		decodePolyId(polyRef, salt, it, ip)
		if (it >= m_maxTiles as UInt) return false
		if (m_tiles[it] salt != salt || m_tiles[it] header == 0) return false
		header: DTMeshHeader* = m_tiles[it] header
		if (ip >= (header polyCount as UInt)) return false
		poly: DTPoly* = &header polys[ip]
		
		// Make sure that the current poly is indeed off-mesh link.
		if ((poly flags & DT_POLY_OFFMESH_CONNECTION) == 0)
			return false
		
		// Figure out which way to hand out the vertices.
		idx0 := 0
		idx1 := 1
		
		// Find link that points to first vertex.
		for (i := poly firstLink; i != DT_NULL_LINK; i = header links[i] next) {
			if (header links[i] egde == 0) {
				if (header links[i] ref != prevRef) {
					idx0 = 1
					idx1 = 0
				}
				break
			}
		}
		vcopy(startPos, &header verts[poly verts[idx0]*3])
		vcopy(endPos, &header verts[poly verts[idx1]*3])
		
		return true
	}
	
	getPolyHeight: func(ref: DTPolyRef, pos, height: Float*) -> Bool {
		salt, it, ip: UInt
		decodePolyId(ref, salt, it, ip)
		if (it >= m_maxTiles as UInt) return false
		if (m_tiles[it] salt != salt || m_tiles[it] header == 0) return false
		header: DTMeshHeader* = m_tiles[it] header
		
		if (ip >= (header polyCount as UInt)) return false
		poly: DTPoly* = &header polys[ip]
		
		if (poly flags & DT_POLY_OFFMESH_CONNECTION) {
			v0: Float* = &header verts[poly verts[0]*3]
			v1: Float* = &header verts[poly verts[1]*3]
			d0: Float = vdist(pos, v0)
			d1: Float = vdist(pos, v1)
			Float u = d0 / (d0+d1)
			if (height)
				*height = v0[1] + (v1[1] - v0[1]) * u
			return true
		} else {
			pd: DTPolyDetail* = &header detailMeshes[ip]
			for (Int j = 0; j < pd triCount; ++j) {
				t: UInt8* = &header detailTris[(pd triBase+j)*4]
				v: Float*[3]
				for (k: Int in 0..3) {
					if (t[k] < poly vertCount)
						v[k] = &header verts[poly verts[t[k]]*3]
					else
						v[k] = &header detailVerts[(pd vertBase+(t[k]-poly vertCount))*3]
				}
				h: Float
				if (closestHeightPointTriangle(pos, v[0], v[1], v[2], h)) {
					if (height)
						*height = h
					return true
				}
			}
		}
		return false
	}
	
	DTPolyRef findNearestPoly(center, extents: Float*, filter: DTQueryFilter*, nearestPt: Float*) {
		// Get nearby polygons from proximity grid.
		polys: DTPolyRef[128]
		Int polyCount = queryPolygons(center, extents, filter, polys, 128)
		
		// Find nearest polygon amongst the nearby polygons.
		nearest: DTPolyRef = 0
		nearestDistanceSqr := FLT_MAX
		for (i: Int in 0..polyCount) {
			ref: DTPolyRef = polys[i]
			closestPtPoly: Float[3]
			if (!closestPointOnPoly(ref, center, closestPtPoly))
				continue
			d: Float = vdistSqr(center, closestPtPoly)
			if (d < nearestDistanceSqr) {
				if (nearestPt)
					vcopy(nearestPt, closestPtPoly)
				nearestDistanceSqr = d
				nearest = ref
			}
		}
		return nearest
	}

	findNearestPolyInTile: func(tile: DTMeshTile*, center, extents: Float*, filter: DTQueryFilter*, nearestPt: Float*) -> DTPolyRef {
		bmin, bmax: Float[3]
		vsub(bmin, center, extents)
		vadd(bmax, center, extents)
		
		// Get nearby polygons from proximity grid.
		polys: DTPolyRef[128]
		polyCount := queryPolygonsInTile(tile, bmin, bmax, filter, polys, 128)
		
		// Find nearest polygon amongst the nearby polygons.
		nearest: DTPolyRef = 0
		nearestDistanceSqr := FLT_MAX
		for (i: Int in 0..polyCount) {
			ref: DTPolyRef = polys[i]
			closestPtPoly: Float[3]
			if (!closestPointOnPolyInTile(tile, decodePolyIdPoly(ref), center, closestPtPoly))
				continue
			d: Float = vdistSqr(center, closestPtPoly)
			if (d < nearestDistanceSqr) {
				if (nearestPt)
					vcopy(nearestPt, closestPtPoly)
				nearestDistanceSqr = d
				nearest = ref
			}
		}
	
		return nearest
	}

	queryPolygonsInTile: func(tile: DTMeshTile*, qmin, qmax: Float*, filter: DTQueryFilter*, DTPolyRef* polys, maxPolys: Int) -> Int {
		header: DTMeshHeader* = tile header
		if (header bvTree) {
			node: DTBVNode* = &header bvTree[0]
			end: DTBVNode* = &header bvTree[header bvNodeCount]
			
			// Calculate quantized box
			bmin, bmax: UShort[3]
			// Clamp query box to world box.
			minx := clamp(qmin[0], header bmin[0], header bmax[0]) - header bmin[0]
			miny := clamp(qmin[1], header bmin[1], header bmax[1]) - header bmin[1]
			minz := clamp(qmin[2], header bmin[2], header bmax[2]) - header bmin[2]
			maxx := clamp(qmax[0], header bmin[0], header bmax[0]) - header bmin[0]
			maxy := clamp(qmax[1], header bmin[1], header bmax[1]) - header bmin[1]
			maxz := clamp(qmax[2], header bmin[2], header bmax[2]) - header bmin[2]
			// Quantize
			bmin[0] = ((header bvQuantFactor * minx) as UShort) & 0xfffe
			bmin[1] = ((header bvQuantFactor * miny) as UShort) & 0xfffe
			bmin[2] = ((header bvQuantFactor * minz) as UShort) & 0xfffe
			bmax[0] = ((header bvQuantFactor * maxx + 1) as UShort) | 1
			bmax[1] = ((header bvQuantFactor * maxy + 1) as UShort) | 1
			bmax[2] = ((header bvQuantFactor * maxz + 1) as UShort) | 1
			
			// Traverse tree
			base: DTPolyRef = getTileId(tile)
			n := 0
			while (node < end) {
				overlap := checkOverlapBox(bmin, bmax, node bmin, node bmax)
				isLeafNode := (node i >= 0)
				
				if (isLeafNode && overlap) {
					if (passFilter(filter, header polys[node i] flags)) {
						if (n < maxPolys)
							polys[n++] = base | (DTPolyRef)node i
					}
				}
			
				if (overlap || isLeafNode)
					node++
				else {
					Int escapeIndex = -node i
					node += escapeIndex
				}
			}
			return n
		} else {
			bmin, bmax: Float[3]
			header: DTMeshHeader* = tile header
			n := 0
			base: DTPolyRef = getTileId(tile)
			for (i: Int in 0..(header polyCount)) {
				// Calc polygon bounds.
				p: DTPoly* = &header polys[i]
				v: Float* = &header verts[p verts[0]*3]
				vcopy(bmin, v)
				vcopy(bmax, v)
				for (j: Int in 1..(p vertCount)) {
					v = &header verts[p verts[j]*3]
					vmin(bmin, v)
					vmax(bmax, v)
				}
				if (overlapBoxes(qmin, qmax, bmin, bmax)) {
					if (passFilter(filter, p flags)) {
						if (n < maxPolys)
							polys[n++] = base | (DTPolyRef)i
					}
				}
			}
			return n
		}
	}
	
	queryPolygons: func(center, extents: Float*, filter: DTQueryFilter*, DTPolyRef* polys, maxPolys: Int) -> Int {
		bmin, bmax: Float[3]
		vsub(bmin, center, extents)
		vadd(bmax, center, extents)
		
		// Find tiles the query touches.
		minx := ((bmin[0]-m_orig[0]) / m_tileWidth) as Int
		maxx := ((bmax[0]-m_orig[0]) / m_tileWidth) as Int
		miny := ((bmin[2]-m_orig[2]) / m_tileHeight) as Int
		maxy := ((bmax[2]-m_orig[2]) / m_tileHeight) as Int
		
		n := 0
		for (y := miny; y <= maxy; ++y) {
			for (x := minx; x <= maxx; ++x) {
				tile: DTMeshTile* = getTileAt(x, y)
				if (!tile) continue
				n += queryPolygonsInTile(tile, bmin, bmax, filter, polys+n, maxPolys-n)
				if (n >= maxPolys) return n
			}
		}
		return n
	}
	
	findPath: func(startRef, endRef: DTPolyRef, startPos, endPos: Float*, filter: DTQueryFilter*, path: DTPolyRef*, maxPathSize: Int) -> Int {
		if (!startRef || !endRef)
			return 0
		
		if (!maxPathSize)
			return 0
		
		if (!getPolyByRef(startRef) || !getPolyByRef(endRef))
			return 0
		
		if (startRef == endRef) {
			path[0] = startRef
			return 1
		}
		
		if (!m_nodePool || !m_openList)
			return 0
		
		m_nodePool clear()
		m_openList clear()
	
		H_SCALE: static const Float = 1.1		// Heuristic scale.
		
		startNode: DTNode* = m_nodePool getNode(startRef)
		startNode pidx = 0
		startNode cost = 0
		startNode total = vdist(startPos, endPos) * H_SCALE
		startNode id = startRef
		startNode flags = DT_NODE_OPEN
		m_openList push(startNode)
	
		lastBestNode: DTNode* = startNode
		lastBestNodeCost := startNode total
		while (!m_openList empty()) {
			bestNode: DTNode* = m_openList pop()
			if (bestNode id == endRef) {
				lastBestNode = bestNode
				break
			}
			
			// Get poly and tile.
			salt, it, ip: UInt
			decodePolyId(bestNode id, salt, it, ip)
			// The API input has been cheked already, skip checking internal data.
			header: DTMeshHeader* = m_tiles[it] header
			poly: DTPoly* = &header polys[ip]
		
			for (i: UInt = poly firstLink; i != DT_NULL_LINK; i = header links[i] next) {
				neighbour: DTPolyRef = header links[i] ref
				if (neighbour) {
					// Skip parent node.
					if (bestNode pidx && m_nodePool getNodeAtIdx(bestNode pidx) id == neighbour)
						continue

					// TODO: Avoid digging the polygon (done in getEdgeMidPoint too).
					if (!passFilter(filter, getPolyFlags(neighbour)))
						continue

					parent: DTNode* = bestNode
					newNode: DTNode
					newNode pidx = m_nodePool getNodeIdx(parent)
					newNode id = neighbour

					// Calculate cost.
					p0, p1: Float[3]
					if (!parent pidx)
						vcopy(p0, startPos)
					else
						getEdgeMidPoint(m_nodePool getNodeAtIdx(parent pidx) id, parent id, p0)
				
					getEdgeMidPoint(parent id, newNode id, p1)
					
					newNode cost = parent cost + vdist(p0, p1)
					// Special case for last node.
					if (newNode id == endRef)
						newNode cost += vdist(p1, endPos)
					
					// Heuristic
					h := vdist(p1, endPos)*H_SCALE
					newNode total = newNode cost + h
					
					actualNode: DTNode* = m_nodePool getNode(newNode id)
					if (!actualNode)
						continue
				
					if (!((actualNode flags & DT_NODE_OPEN) && newNode total > actualNode total) &&
						!((actualNode flags & DT_NODE_CLOSED) && newNode total > actualNode total)) {
						actualNode flags &= ~DT_NODE_CLOSED
						actualNode pidx = newNode pidx
						actualNode cost = newNode cost
						actualNode total = newNode total
						
						if (h < lastBestNodeCost) {
							lastBestNodeCost = h
							lastBestNode = actualNode
						}
						if (actualNode flags & DT_NODE_OPEN) {
							m_openList modify(actualNode)
						} else {
							actualNode flags |= DT_NODE_OPEN
							m_openList push(actualNode)
						}
					}
				}
			}
			bestNode flags |= DT_NODE_CLOSED
		}
		
		// Reverse the path.
		prev: DTNode* = 0
		node: DTNode* = lastBestNode
		while (true) { //do
			next: DTNode* = m_nodePool getNodeAtIdx(node pidx)
			node pidx = m_nodePool getNodeIdx(prev)
			prev = node
			node = next
			if (!node) break
		} //while (node)
		
		// Store path
		node = prev
		n := 0
		while (true) { //do
			path[n++] = node id
			node = m_nodePool getNodeAtIdx(node pidx)
			if (!node || n >= maxPathSize) break
		} //while (node && n < maxPathSize)
		return n
	}

	findStraightPath: func(startPos, endPos: Float*, path: DTPolyRef*, pathSize: Int, straightPath: Float*,
							straightPathFlags: UInt8*, straightPathRefs: DTPolyRef*, maxStraightPathSize: Int) -> Int {
		if (!maxStraightPathSize)
			return 0
		if (!path[0])
			return 0
		
		straightPathSize := 0
		
		// TODO: Should this be callers responsibility?
		closestStartPos: Float[3]
		if (!closestPointOnPolyBoundary(path[0], startPos, closestStartPos))
			return 0
		
		// Add start point.
		vcopy(&straightPath[straightPathSize*3], closestStartPos)
		if (straightPathFlags)
			straightPathFlags[straightPathSize] = DT_STRAIGHTPATH_START
		if (straightPathRefs)
			straightPathRefs[straightPathSize] = path[0]
		straightPathSize++
		if (straightPathSize >= maxStraightPathSize)
			return straightPathSize
	
		closestEndPos: Float[3]
		if (!closestPointOnPolyBoundary(path[pathSize-1], endPos, closestEndPos))
			return 0
		
		if (pathSize > 1) {
			Float portalApex[3], portalLeft[3], portalRight[3]
			vcopy(portalApex, closestStartPos)
			vcopy(portalLeft, portalApex)
			vcopy(portalRight, portalApex)
			apexIndex := 0
			leftIndex := 0
			rightIndex := 0
			
			leftPolyFlags: UShort = 0
			rightPolyFlags: UShort = 0
			
			leftPolyRef: DTPolyRef = path[0]
			rightPolyRef: DTPolyRef = path[0]
			
			for (i: Int in 0..pathSize) {
				left, right: Float[3]
				fromFlags, toFlags: UShort
				
				if (i+1 < pathSize) {
					// Next portal.
					if (!getPortalPoints(path[i], path[i+1], left, right, fromFlags, toFlags)) {
						if (!closestPointOnPolyBoundary(path[i], endPos, closestEndPos))
							return 0
						
						vcopy(&straightPath[straightPathSize*3], closestEndPos)
						if (straightPathFlags)
							straightPathFlags[straightPathSize] = 0
						if (straightPathRefs)
							straightPathRefs[straightPathSize] = path[i]
						straightPathSize++
					
						return straightPathSize
					}
				} else {
					// End of the path.
					vcopy(left, closestEndPos)
					vcopy(right, closestEndPos)

					fromFlags = toFlags = 0
				}
			
				// Right vertex.
				if (vequal(portalApex, portalRight)) {
					vcopy(portalRight, right)
					rightPolyRef = (i+1 < pathSize) ? path[i+1] : 0
					rightPolyFlags = toFlags
					rightIndex = i
				} else {
					if (triArea2D(portalApex, portalRight, right) <= 0.0) {
						if (triArea2D(portalApex, portalLeft, right) > 0.0) {
							vcopy(portalRight, right)
							rightPolyRef = (i+1 < pathSize) ? path[i+1] : 0
							rightPolyFlags = toFlags
							rightIndex = i
						} else {
							vcopy(portalApex, portalLeft)
							apexIndex = leftIndex
						
							flags: UInt8 = (leftPolyFlags & DT_POLY_OFFMESH_CONNECTION) ? DT_STRAIGHTPATH_OFFMESH_CONNECTION : 0
							ref: DTPolyRef = leftPolyRef
							
							if (!vequal(&straightPath[(straightPathSize-1)*3], portalApex)) {
								vcopy(&straightPath[straightPathSize*3], portalApex)
								if (straightPathFlags)
									straightPathFlags[straightPathSize] = flags
								if (straightPathRefs)
									straightPathRefs[straightPathSize] = ref
								
								straightPathSize++
								if (straightPathSize >= maxStraightPathSize)
									return straightPathSize
							} else {
								// The vertices are equal, update flags and poly.
								if (straightPathFlags)
									straightPathFlags[straightPathSize-1] = flags
								if (straightPathRefs)
									straightPathRefs[straightPathSize-1] = ref
							}
						
							vcopy(portalLeft, portalApex)
							vcopy(portalRight, portalApex)
							leftIndex = apexIndex
							rightIndex = apexIndex
							// Restart
							i = apexIndex
							continue
						}
					}
				}
				
				// Left vertex.
				if (vequal(portalApex, portalLeft)) {
					vcopy(portalLeft, left)
					leftPolyRef = (i+1 < pathSize) ? path[i+1] : 0
					leftPolyFlags = toFlags
					leftIndex = i
				} else {
					if (triArea2D(portalApex, portalLeft, left) >= 0.0) {
						if (triArea2D(portalApex, portalRight, left) < 0.0) {
							vcopy(portalLeft, left)
							leftPolyRef = (i+1 < pathSize) ? path[i+1] : 0
							leftPolyFlags = toFlags
							leftIndex = i
						} else {
							vcopy(portalApex, portalRight)
							apexIndex = rightIndex
							
							flags: UInt8 = (rightPolyFlags & DT_POLY_OFFMESH_CONNECTION) ? DT_STRAIGHTPATH_OFFMESH_CONNECTION : 0
							ref: DTPolyRef = rightPolyRef
						
							if (!vequal(&straightPath[(straightPathSize-1)*3], portalApex)) {
								vcopy(&straightPath[straightPathSize*3], portalApex)
								if (straightPathFlags)
									straightPathFlags[straightPathSize] = flags
								if (straightPathRefs)
									straightPathRefs[straightPathSize] = ref
							
								straightPathSize++
								if (straightPathSize >= maxStraightPathSize)
									return straightPathSize
							} else {
								// The vertices are equal, update flags and poly.
								if (straightPathFlags)
									straightPathFlags[straightPathSize-1] = flags
								if (straightPathRefs)
									straightPathRefs[straightPathSize-1] = ref
							}
							
							vcopy(portalLeft, portalApex)
							vcopy(portalRight, portalApex)
							leftIndex = apexIndex
							rightIndex = apexIndex
							
							// Restart
							i = apexIndex
						
							continue
						}
					}
				}
			}
		}
		
		// Add end point.
		vcopy(&straightPath[straightPathSize*3], closestEndPos)
		if (straightPathFlags)
			straightPathFlags[straightPathSize] = DT_STRAIGHTPATH_END
		if (straightPathRefs)
			straightPathRefs[straightPathSize] = 0
		
		straightPathSize++
		return straightPathSize
	}

	// Moves towards end position a long the path corridor.
	// Returns: Index to the result path polygon.
	moveAlongPathCorridor: func(startPos, endPos, resultPos: Float*, path: DTPolyRef*, pathSize: Int) -> Int {
		if (!pathSize)
			return 0
		
		verts: Float[DT_VERTS_PER_POLYGON*3]	
		edged: Float[DT_VERTS_PER_POLYGON]
		edget: Float[DT_VERTS_PER_POLYGON]
		n := 0
		
		SLOP: static const Float = 0.01
		vcopy(resultPos, startPos)
		while (n < pathSize) {
			// Get current polygon and poly vertices.
			salt, it, ip: UInt
			decodePolyId(path[n], salt, it, ip)
			if (it >= m_maxTiles as UInt) return n
			if (m_tiles[it] salt != salt || m_tiles[it] header == 0) return n
			if (ip >= (m_tiles[it] header polyCount as UInt)) return n
			header: DTMeshHeader* = m_tiles[it] header
			poly: DTPoly* = &header polys[ip]
			
			// In case of Off-Mesh link, just snap to the end location and advance over it.
			if (poly flags & DT_POLY_OFFMESH_CONNECTION) {
				if (n+1 < pathSize) {
					left, right: Float[3]
					fromFlags, toFlags: UShort
					if (!getPortalPoints(path[n], path[n+1], left, right, fromFlags, toFlags))
						return n
					vcopy(resultPos, endPos)
				}
				return n+1
			}
			
			// Collect vertices.
			nv := 0
			for (i: Int in 0..(poly vertCount as Int)) {
				vcopy(&verts[nv*3], &header verts[poly verts[i]*3])
				nv++
			}
			
			inside := distancePtPolyEdgesSqr(endPos, verts, nv, edged, edget)
			if (inside) {
				// The end point is inside the current polygon.
				vcopy(resultPos, endPos)
				return n
			}

			// Constraint the point on the polygon boundary.
			// This results sliding movement.
			dmin := FLT_MAX
			imin := -1
			for (i: Int in 0..nv) {
				if (edged[i] < dmin) {
					dmin = edged[i]
					imin = i
				}
			}
			va: Float* = &verts[imin*3]
			vb: Float* = &verts[((imin+1)%nv)*3]
			vlerp(resultPos, va, vb, edget[imin])
		
			// Check to see if the point is on the portal edge to the next polygon.
			if (n+1 >= pathSize)
				return n
			left, right: Float[3]
			fromFlags, toFlags: UShort
			if (!getPortalPoints(path[n], path[n+1], left, right, fromFlags, toFlags))
				return n
			// If the clamped point is close to the next portal edge, advance to next poly.
			t: Float
			d := distancePtSegSqr2D(resultPos, left, right, t)
			if (d > SLOP*SLOP)
				return n
			// Advance to next polygon.
			n++
		}
		return n
	}

	// Returns portal points between two polygons.
	getPortalPoints: func(frompr, to: DTPolyRef, left, right: Float*, fromFlags, toFlags: UShort&) -> Bool {
		salt, it, ip: UInt
		decodePolyId(frompr, salt, it, ip)
		if (it >= m_maxTiles as UInt) return false
		if (m_tiles[it] salt != salt || m_tiles[it] header == 0) return false
		if (ip >= (m_tiles[it] header polyCount as UInt)) return false
		fromHeader: DTMeshHeader* = m_tiles[it] header
		fromPoly: DTPoly* = &fromHeader polys[ip]
		fromFlags = fromPoly flags
		
		for (i: UInt = fromPoly firstLink; i != DT_NULL_LINK; i = fromHeader links[i] next) {
			link: DTLink* = &fromHeader links[i]
			if (link ref != to)
				continue

			decodePolyId(to, salt, it, ip)
			if (it >= m_maxTiles as UInt) return false
			if (m_tiles[it] salt != salt || m_tiles[it] header == 0) return false
			if (ip >= (m_tiles[it] header polyCount as UInt)) return false
			toHeader: DTMeshHeader* = m_tiles[it] header
			toPoly: DTPoly* = &toHeader polys[ip]
			toFlags = toPoly flags
			
			if (fromPoly flags & DT_POLY_OFFMESH_CONNECTION) {
				// Find link that points to first vertex.
				for (i: UInt = fromPoly firstLink; i != DT_NULL_LINK; i = fromHeader links[i] next) {
					if (fromHeader links[i] ref == to) {
						v := fromHeader links[i] egde
						vcopy(left, &fromHeader verts[fromPoly verts[v]*3])
						vcopy(right, &fromHeader verts[fromPoly verts[v]*3])
						return true
					}
				}
				return false
			}

			if (toPoly flags & DT_POLY_OFFMESH_CONNECTION) {
				for (i: UInt = toPoly firstLink; i != DT_NULL_LINK; i = toHeader links[i] next) {
					if (toHeader links[i] ref == frompr) {
						v := toHeader links[i] egde
						vcopy(left, &toHeader verts[toPoly verts[v]*3])
						vcopy(right, &toHeader verts[toPoly verts[v]*3])
						return true
					}
				}
				return false
			}
		
			// Find portal vertices.
			v0 := fromPoly verts[link edge]
			v1 := fromPoly verts[(link edge+1) % (fromPoly vertCount as Int)]
			vcopy(left, &fromHeader verts[v0*3])
			vcopy(right, &fromHeader verts[v1*3])
			// If the link is at tile boundary, clamp the vertices to
			// the link width.
			if (link side == 0 || link side == 4) {
				// Unpack portal limits.
				smin := min(left[2], right[2])
				smax := max(left[2], right[2])
				s: Float = (smax-smin) / 255.0f
				lmin: Float = smin + link bmin*s
				lmax: Float = smin + link bmax*s
				left[2] = max(left[2], lmin)
				left[2] = min(left[2], lmax)
				right[2] = max(right[2], lmin)
				right[2] = min(right[2], lmax)
			} else if (link side == 2 || link side == 6) {
				// Unpack portal limits.
				smin := min(left[0], right[0])
				smax := max(left[0], right[0])
				s: Float = (smax-smin) / 255.0f
				lmin: Float = smin + link bmin*s
				lmax: Float = smin + link bmax*s
				left[0] = max(left[0], lmin)
				left[0] = min(left[0], lmax)
				right[0] = max(right[0], lmin)
				right[0] = min(right[0], lmax)
			}
			return true
		}
		return false
	}
	
	// Returns edge mid point between two polygons.
	getEdgeMidPoint: func(frompr, to: DTPolyRef, Float* mid) -> Bool {
		left, right: Float[3]
		fromFlags, toFlags: UShort
		if (!getPortalPoints(frompr, to, left, right, fromFlags, toFlags)) return false
		mid[0] = (left[0]+right[0])*0.5f
		mid[1] = (left[1]+right[1])*0.5f
		mid[2] = (left[2]+right[2])*0.5f
		return true
	}

	getPolyFlags: func(ref: DTPolyRef) -> UShort {
		salt, it, ip: UInt
		decodePolyId(ref, salt, it, ip)
		if (it >= m_maxTiles as UInt) return 0
		if (m_tiles[it] salt != salt || m_tiles[it] header == 0) return 0
		if (ip >= (m_tiles[it] header polyCount as UInt)) return 0
		header: DTMeshHeader* = m_tiles[it] header
		poly: DTPoly* = &header polys[ip]
		return poly flags
	}
	
	raycast: func(centerRef: DTPolyRef, startPos, endPos: Float*, filter: DTQueryFilter*,
					t: Float&, hitNormal: Float*, path: DTPolyRef*, pathSize: Int) -> Int {
		t = 0
	
		if (!centerRef || !getPolyByRef(centerRef))
			return 0
	
		DTPolyRef curRef = centerRef
		verts: Float[DT_VERTS_PER_POLYGON*3];	
		n := 0
	
		hitNormal[0] = 0
		hitNormal[1] = 0
		hitNormal[2] = 0
	
		while (curRef) {
			// Cast ray against current polygon.
		
			// The API input has been cheked already, skip checking internal data.
			UInt it = decodePolyIdTile(curRef)
			UInt ip = decodePolyIdPoly(curRef)
			header: DTMeshHeader* = m_tiles[it] header
			poly: DTPoly* = &header polys[ip]

			// Collect vertices.
			nv := 0
			for (i: Int in 0..(poly vertCount as Int)) {
				vcopy(&verts[nv*3], &header verts[poly verts[i]*3])
				nv++
			}		
		
			Float tmin, tmax
			Int segMin, segMax
			if (!intersectSegmentPoly2D(startPos, endPos, verts, nv, tmin, tmax, segMin, segMax)) {
				// Could not hit the polygon, keep the old t and report hit.
				return n
			}
			// Keep track of furthest t so far.
			if (tmax > t)
				t = tmax

			if (n < pathSize)
				path[n++] = curRef
		
			// Follow neighbours.
			nextRef: DTPolyRef = 0
			
			for (i: UInt = poly firstLink; i != DT_NULL_LINK; i = header links[i] next) {
				link: DTLink* = &header links[i]
				if ((Int)link edge == segMax) {
					// If the link is internal, just return the ref.
					if (link side == 0xff) {
						nextRef = link ref
						break
					}
					
					// If the link is at tile boundary, 
					v0 := poly verts[link edge]
					v1 := poly verts[(link edge+1) % poly vertCount]
					left: Float* = &header verts[v0*3]
					right: Float* = &header verts[v1*3]
				
					// Check that the intersection lies inside the link portal.
					if (link side == 0 || link side == 4) {
						// Calculate link size.
						smin := min(left[2], right[2])
						smax := max(left[2], right[2])
						s: Float = (smax-smin) / 255.0
						lmin: Float = smin + link bmin*s
						lmax: Float = smin + link bmax*s
						// Find Z intersection.
						z: Float = startPos[2] + (endPos[2]-startPos[2])*tmax
						if (z >= lmin && z <= lmax) {
							nextRef = link ref
							break
						}
					} else if (link side == 2 || link side == 6) {
						// Calculate link size.
						smin := min(left[0], right[0])
						smax := max(left[0], right[0])
						s: Float = (smax-smin) / 255.0
						lmin: Float = smin + link bmin*s
						lmax: Float = smin + link bmax*s
						// Find X intersection.
						x: Float = startPos[0] + (endPos[0]-startPos[0])*tmax
						if (x >= lmin && x <= lmax) {
							nextRef = link ref
							break
						}
					}
				}
			}
		
			if (!nextRef || !passFilter(filter, getPolyFlags(nextRef))) {
				// No neighbour, we hit a wall.

				// Calculate hit normal.
				a := segMax
				b := segMax+1 < nv ? segMax+1 : 0
				va: Float* = &verts[a*3]
				vb: Float* = &verts[b*3]
				d: Floatx = vb[0] - va[0]
				d: Floatz = vb[2] - va[2]
				hitNormal[0] = dz
				hitNormal[1] = 0
				hitNormal[2] = -dx
				vnormalize(hitNormal)
				return n
			}
		
			// No hit, advance to neighbour polygon.
			curRef = nextRef
		}
		return n
	}
	
	findPolysAround: func(centerRef: DTPolyRef, centerPos: Float*, radius: Float, filter: DTQueryFilter*, 
							resultRef, resultParent: DTPolyRef*, resultCost: Float*, maxResult: Int) -> Int {
		if (!centerRef) return 0
		if (!getPolyByRef(centerRef)) return 0
		if (!m_nodePool || !m_openList) return 0
		
		m_nodePool clear()
		m_openList clear()
		
		startNode: DTNode* = m_nodePool getNode(centerRef)
		startNode pidx = 0
		startNode cost = 0
		startNode total = 0
		startNode id = centerRef
		startNode flags = DT_NODE_OPEN
		m_openList push(startNode)
		
		n := 0
		if (n < maxResult) {
			if (resultRef)
				resultRef[n] = startNode id
			if (resultParent)
				resultParent[n] = 0
			if (resultCost)
				resultCost[n] = 0
			++n
		}
		
		radiusSqr: Float = sqr(radius)
		while (!m_openList empty()) {
			bestNode: DTNode* = m_openList pop()

			// Get poly and tile.
			// The API input has been cheked already, skip checking internal data.
			UInt it = decodePolyIdTile(bestNode id)
			UInt ip = decodePolyIdPoly(bestNode id)
			header: DTMeshHeader* = m_tiles[it] header
			poly: DTPoly* = &header polys[ip]
		
			for (i: UInt = poly firstLink; i != DT_NULL_LINK; i = header links[i] next) {
				link: DTLink* = &header links[i]
				neighbour: DTPolyRef = link ref
				if (neighbour) {
					// Skip parent node.
					if (bestNode pidx && m_nodePool getNodeAtIdx(bestNode pidx) id == neighbour)
						continue
					
					// Calc distance to the edge.
					va: Float* = &header verts[poly verts[link edge]*3]
					vb: Float* = &header verts[poly verts[(link edge+1) % poly vertCount]*3]
					tseg: Float
					d: FloatistSqr = distancePtSegSqr2D(centerPos, va, vb, tseg)
					
					// If the circle is not touching the next polygon, skip it.
					if (distSqr > radiusSqr)
						continue
					if (!passFilter(filter, getPolyFlags(neighbour)))
						continue
					
					parent: DTNode* = bestNode
					newNode: DTNode
					newNode pidx = m_nodePool getNodeIdx(parent)
					newNode id = neighbour

					// Cost
					p0, p1: Float[3]
					if (!parent pidx)
						vcopy(p0, centerPos)
					else
						getEdgeMidPoint(m_nodePool getNodeAtIdx(parent pidx) id, parent id, p0)
					getEdgeMidPoint(parent id, newNode id, p1)
					newNode total = parent total + vdist(p0, p1)
				
					actualNode: DTNode* = m_nodePool getNode(newNode id)
					if (!actualNode)
						continue
					
					if (!((actualNode flags & DT_NODE_OPEN) && newNode total > actualNode total) &&
						!((actualNode flags & DT_NODE_CLOSED) && newNode total > actualNode total)) {
						actualNode flags &= ~DT_NODE_CLOSED
						actualNode pidx = newNode pidx
						actualNode total = newNode total
					
						if (actualNode flags & DT_NODE_OPEN) {
							m_openList modify(actualNode)
						} else {
							if (n < maxResult) {
								if (resultRef)
									resultRef[n] = actualNode id
								if (resultParent)
									resultParent[n] = m_nodePool getNodeAtIdx(actualNode pidx) id
								if (resultCost)
									resultCost[n] = actualNode total
								++n
							}
							actualNode flags = DT_NODE_OPEN
							m_openList push(actualNode)
						}
					}
				}
			}
		}
		return n
	}

	findDistanceToWall: func(centerRef: DTPolyRef, centerPos: Float*, maxRadius: Float, filter: DTQueryFilter*, hitPos, hitNormal: Float*) -> Float {
		if (!centerRef) return 0
		if (!getPolyByRef(centerRef)) return 0
		if (!m_nodePool || !m_openList) return 0
		
		m_nodePool clear()
		m_openList clear()
	
		startNode: DTNode* = m_nodePool getNode(centerRef)
		startNode pidx = 0
		startNode cost = 0
		startNode total = 0
		startNode id = centerRef
		startNode flags = DT_NODE_OPEN
		m_openList push(startNode)
	
		radiusSqr: Float = sqr(maxRadius)
	
		while (!m_openList empty()) {
			bestNode: DTNode* = m_openList pop()
		
			// Get poly and tile.
			// The API input has been cheked already, skip checking internal data.
			it := decodePolyIdTile(bestNode id)
			ip := decodePolyIdPoly(bestNode id)
			header: DTMeshHeader* = m_tiles[it] header
			poly: DTPoly* = &header polys[ip]
		
			// Hit test walls.
			j := (poly vertCount as Int) - 1
			for (i: Int in 0..(poly vertCount as Int)) {
				// Skip non-solid edges.
				if (poly neis[j] & DT_EXT_LINK) {
					// Tile border.
					solid := true
					for (k: UInt = poly firstLink; k != DT_NULL_LINK; k = header links[k] next) {
						link: DTLink* = &header links[k]
						if (link edge == j && link ref != 0 && passFilter(filter, getPolyFlags(link ref))) {
							solid = false
							break
						}
					}
					if (!solid) continue
				} else if (poly neis[j] && passFilter(filter, header polys[poly neis[j]] flags)) {
					// Internal edge
					continue
				}
			
				// Calc distance to the edge.
				vj: Float* = &header verts[poly verts[j]*3]
				vi: Float* = &header verts[poly verts[i]*3]
				tseg: Float
				d: FloatistSqr = distancePtSegSqr2D(centerPos, vj, vi, tseg)
			
				// Edge is too far, skip.
				if (distSqr > radiusSqr)
					continue
			
				// Hit wall, update radius.
				radiusSqr = distSqr
				// Calculate hit pos.
				hitPos[0] = vj[0] + (vi[0] - vj[0])*tseg
				hitPos[1] = vj[1] + (vi[1] - vj[1])*tseg
				hitPos[2] = vj[2] + (vi[2] - vj[2])*tseg
			}
			
			for (i: UInt = poly firstLink; i != DT_NULL_LINK; i = header links[i] next) {
				link: DTLink* = &header links[i]
				neighbour: DTPolyRef = link ref
				if (neighbour) {
					// Skip parent node.
					if (bestNode pidx && m_nodePool getNodeAtIdx(bestNode pidx) id == neighbour)
						continue
					
					// Calc distance to the edge.
					va: Float* = &header verts[poly verts[link edge]*3]
					vb: Float* = &header verts[poly verts[(link edge+1) % poly vertCount]*3]
					tseg: Float
					d: FloatistSqr = distancePtSegSqr2D(centerPos, va, vb, tseg)
					
					// If the circle is not touching the next polygon, skip it.
					if (distSqr > radiusSqr)
						continue
					if (!passFilter(filter, getPolyFlags(neighbour)))
						continue
					
					parent: DTNode* = bestNode
					newNode: DTNode
					newNode pidx = m_nodePool getNodeIdx(parent)
					newNode id = neighbour
					
					p0, p1: Float[3]
					if (!parent pidx)
						vcopy(p0, centerPos)
					else
						getEdgeMidPoint(m_nodePool getNodeAtIdx(parent pidx) id, parent id, p0)
					getEdgeMidPoint(parent id, newNode id, p1)
					newNode total = parent total + vdist(p0, p1)
					
					actualNode: DTNode* = m_nodePool getNode(newNode id)
					if (!actualNode)
						continue
					
					if (!((actualNode flags & DT_NODE_OPEN) && newNode total > actualNode total) &&
						!((actualNode flags & DT_NODE_CLOSED) && newNode total > actualNode total)) {
						actualNode flags &= ~DT_NODE_CLOSED
						actualNode pidx = newNode pidx
						actualNode total = newNode total
						
						if (actualNode flags & DT_NODE_OPEN) {
							m_openList modify(actualNode)
						} else {
							actualNode flags = DT_NODE_OPEN
							m_openList push(actualNode)
						}
					}
				}
			}
		}
		
		// Calc hit normal.
		vsub(hitNormal, centerPos, hitPos)
		vnormalize(hitNormal)
		return sqrtf(radiusSqr)
	}
	
	getPolyByRef: func(ref: DTPolyRef) -> DTPoly* {
		salt, it, ip: UInt
		decodePolyId(ref, salt, it, ip)
		if (it >= m_maxTiles as UInt) return 0
		if (m_tiles[it] salt != salt || m_tiles[it] header == 0) return 0
		if (ip >= (m_tiles[it] header polyCount as UInt)) return 0
		return &m_tiles[it] header polys[ip]
	}
	
	getPolyVertsByRef: func(ref: DTPolyRef) -> Float* {
		salt, it, ip: UInt
		decodePolyId(ref, salt, it, ip)
		if (it >= m_maxTiles as UInt) return 0
		if (m_tiles[it] salt != salt || m_tiles[it] header == 0) return 0
		if (ip >= (m_tiles[it] header polyCount as UInt)) return 0
		return m_tiles[it] header verts
	}
	
	getPolyLinksByRef: func(ref: DTPolyRef) -> DTLink* {
		salt, it, ip: UInt
		decodePolyId(ref, salt, it, ip)
		if (it >= m_maxTiles as UInt) return 0
		if (m_tiles[it] salt != salt || m_tiles[it] header == 0) return 0
		if (ip >= (m_tiles[it] header polyCount as UInt)) return 0
		return m_tiles[it] header links
	}
	
	isInClosedList: func(ref: DTPolyRef) -> Bool {
		if (!m_nodePool) return false
		node: DTNode* = m_nodePool findNode(ref)
		return node && node flags & DT_NODE_CLOSED
	}
	
	encodePolyId: inline func(salt, it, ip: UInt) -> DTPolyRef {
		return (salt << (m_polyBits+m_tileBits)) | ((it+1) << m_polyBits) | ip
	}
	
	decodePolyId: inline func(ref: DTPolyRef, salt, it, ip: UInt&) {
		salt = (ref >> (m_polyBits+m_tileBits)) & ((1<<m_saltBits)-1)
		it = ((ref >> m_polyBits) - 1) & ((1<<m_tileBits)-1)
		ip = ref & ((1<<m_polyBits)-1)
	}
	
	decodePolyIdTile inline func(ref: DTPolyRef) -> UInt {
		return ((ref >> m_polyBits) - 1) & ((1<<m_tileBits)-1)
	}
	
	decodePolyIdPoly: inline func(ref: DTPolyRef) -> UInt {
		return ref & ((1<<m_polyBits)-1)
	}
	
// *** private
	m_orig: Float[3]					// Origin of the tile (0,0)
	m_tileWidth, m_tileHeight: Float	// Dimensions of each tile.
	m_maxTiles: Int						// Max number of tiles.
	m_tileLutSize: Int					// Tile hash lookup size (must be pot).
	m_tileLutMask: Int					// Tile hash lookup mask.
	
	m_posLookup: DTMeshTile**			// Tile hash lookup.
	m_nextFree: DTMeshTile*				// Freelist of tiles.
	m_tiles: DTMeshTile*				// List of tiles.
	
	m_saltBits: UInt			// Number of salt bits in the tile ID.
	m_tileBits: UInt			// Number of tile bits in the tile ID.
	m_polyBits: UInt			// Number of poly bits in the tile ID.
	
	m_nodePool: DTNodePool*		// Pointer to node pool.
	m_openList: DTNodeQueue*		// Pointer to open list queue.
}

