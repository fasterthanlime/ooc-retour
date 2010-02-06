use retour
import Retour/[Retour, Timer]

RCEdge: class {
	vert: UShort[2]
	polyEdge: UShort[2]
	poly: UShort[2]
}

buildMeshAdjacency: static func(polys: UShort*, npolys, nverts, vertsPerPoly: Int) -> Bool {
	// Based on code by Eric Lengyel from:
	// http://www.terathon.com/code/edges.php
	
	maxEdgeCount := npolys * vertsPerPoly
	firstEdge := rcAllocArray(UShort, nverts + maxEdgeCount)
	if (!firstEdge)
		return false
	nextEdge := firstEdge + nverts
	edgeCount := 0
	
	edges := rcAllocArray(RCEdge, maxEdgeCount)
	if (!edges)
		return false
	
	for (i: Int in 0..nverts)
		firstEdge[i] = RC_MESH_NULL_IDX
	
	for (i: Int in 0..npolys) {
		t: UShort* = polys[i*vertsPerPoly*2]&
		for (j: Int in 0..vertsPerPoly) {
			v0: UShort = t[j]
			v1: UShort = (j+1 >= vertsPerPoly || t[j+1] == RC_MESH_NULL_IDX) ? t[0] : t[j+1]
			if (v0 < v1) {
				edge: RCEdge& = edges[edgeCount]
				edge vert[0] = v0
				edge vert[1] = v1
				edge poly[0] = i as UShort
				edge polyEdge[0] = j as UShort
				edge poly[1] = i as UShort
				edge polyEdge[1] = 0
				// Insert edge
				nextEdge[edgeCount] = firstEdge[v0]
				firstEdge[v0] = edgeCount
				edgeCount+=1
			}
		}
	}
	
	for (i: Int in 0..npolys) {
		t: UShort* = polys[i*vertsPerPoly*2]&
		for (j: Int in 0..vertsPerPoly) {
			v0:UShort = t[j]
			v1: UShort = (j+1 >= vertsPerPoly || t[j+1] == RC_MESH_NULL_IDX) ? t[0] : t[j+1]
			if (v0 > v1) {
				for (e: UShort = firstEdge[v1]; e != RC_MESH_NULL_IDX; e = nextEdge[e]) {
					edge: RCEdge& = edges[e]
					if (edge vert[1] == v0 && edge poly[0] == edge poly[1]) {
						edge poly[1] = i as UShort
						edge polyEdge[1] = j as UShort
						break
					}
				}
			}
		}
	}
	
	// Store adjacency
	for (i: Int in 0..edgeCount) {
		e: RCEdge = edges[i]
		if (e poly[0] != e poly[1]) {
			p0: UShort* = polys[e poly[0]*vertsPerPoly*2]&
			p1: UShort* = polys[e poly[1]*vertsPerPoly*2]&
			p0[vertsPerPoly + e polyEdge[0]] = e poly[1]
			p1[vertsPerPoly + e polyEdge[1]] = e poly[0]
		}
	}
	
	delete [] firstEdge
	delete [] edges
	
	return true
}

VERTEX_BUCKET_COUNT: static const Int = (1<<12)

computeVertexHash: inline func(x, y, z: Int) -> Int {
	h1: UInt = 0x8da6b343	// Large multiplicative constants
	h2: UInt = 0xd8163841	// here arbitrarily chosen primes
	h3: UInt = 0xcb1ab31f
	n: UInt = h1 * x + h2 * y + h3 * z
	return (n & (VERTEX_BUCKET_COUNT-1)) as Int
}

addVertex: static func(x, y, z: UShort, verts: UShort*, firstVert, nextVert: Int*, nv: Int@) -> Int {
	bucket := computeVertexHash(x, 0, z)
	i := firstVert[bucket]
	
	while (i != -1) {
		v: UShort* = verts[i*3]&
		if (v[0] == x && (rcAbs(v[1] - y) <= 2) && v[2] == z)
			return i
		i = nextVert[i] // next
	}
	
	// Could not find, create new.
	i = nv; nv+=1
	v: UShort* = verts[i*3]&
	v[0] = x
	v[1] = y
	v[2] = z
	nextVert[i] = firstVert[bucket]
	firstVert[bucket] = i
	
	return i
}

prev: inline func(i, n: Int) -> Int {
	return i-1 >= 0 ? i-1 : n-1
}

next: inline func(i, n: Int) -> Int {
	return i+1 < n ? i+1 : 0
}

area2: inline func(a, b, c: Int*) -> Int{
	return (b[0] - a[0]) * (c[2] - a[2]) - (c[0] - a[0]) * (b[2] - a[2])
}

//	Exclusive or: true iff exactly one argument is true.
//	The arguments are negated to ensure that they are 0/1
//	values.  Then the bitwise Xor operator may apply.
//	(This idea is due to Michael Baldwin.)
xorb: inline func(x, y: Bool) -> Bool {
	return !x ^ !y
}

// Returns true iff c is strictly to the left of the directed
// line through a to b.
left: inline func(a, b, c: Int*) -> Bool {
	return area2(a, b, c) < 0
}

leftOn: inline func(a, b, c: Int*) -> Bool {
	return area2(a, b, c) <= 0
}

collinear: inline func(a, b, c: Int*) -> Bool {
	return area2(a, b, c) == 0
}

//	Returns true iff ab properly intersects cd: they share
//	a point interior to both segments.  The properness of the
//	intersection is ensured by using strict leftness.
intersectProp: func(a, b, c, d: Int*) -> Bool {
	// Eliminate improper cases.
	if (collinear(a,b,c) || collinear(a,b,d) ||
		collinear(c,d,a) || collinear(c,d,b))
		return false
	
	return xorb(left(a,b,c), left(a,b,d)) && xorb(left(c,d,a), left(c,d,b))
}

// Returns T iff (a,b,c) are collinear and point c lies 
// on the closed segement ab.
between: static func(a, b, c: Int*) -> Bool {
	if (!collinear(a, b, c))
		return false
	// If ab not vertical, check betweenness on x; else on y.
	if (a[0] != b[0])
		return	((a[0] <= c[0]) && (c[0] <= b[0])) || ((a[0] >= c[0]) && (c[0] >= b[0]))
	else
		return	((a[2] <= c[2]) && (c[2] <= b[2])) || ((a[2] >= c[2]) && (c[2] >= b[2]))
}

// Returns true iff segments ab and cd intersect, properly or improperly.
intersect: static func(a, b, c, d: Int*) -> Bool {
	if (intersectProp(a, b, c, d))
		return true
	else if (between(a, b, c) || between(a, b, d) ||
			 between(c, d, a) || between(c, d, b))
		return true
	else
		return false
}

vequal: static func(a, b: Int*) -> Bool {
	return a[0] == b[0] && a[2] == b[2]
}

// Returns T iff (v_i, v_j) is a proper internal *or* external
// diagonal of P, *ignoring edges incident to v_i and v_j*.
diagonalie: static func(i, j, n: Int, verts, indices: Int*) -> Bool {
	d0: Int* = verts[(indices[i] & 0x0fffffff) * 4]&
	d1: Int* = verts[(indices[j] & 0x0fffffff) * 4]&
	
	// For each edge (k,k+1) of P
	for (k: Int in 0..n) {
		k1 := next(k, n)
		// Skip edges incident to i or j
		if (!((k == i) || (k1 == i) || (k == j) || (k1 == j))) {
			p0: Int* = verts[(indices[k] & 0x0fffffff) * 4]&
			p1: Int* = verts[(indices[k1] & 0x0fffffff) * 4]&
			
			if (vequal(d0, p0) || vequal(d1, p0) || vequal(d0, p1) || vequal(d1, p1))
				continue
			if (intersect(d0, d1, p0, p1))
				return false
		}
	}
	return true
}

// Returns true iff the diagonal (i,j) is strictly internal to the 
// polygon P in the neighborhood of the i endpoint.
inCone: static func(i, j, n: Int, verts, indices: Int*) -> Bool {
	pi: Int* = verts[(indices[i] & 0x0fffffff) * 4]&
	pj: Int* = verts[(indices[j] & 0x0fffffff) * 4]&
	pi1: Int* = verts[(indices[next(i, n)] & 0x0fffffff) * 4]&
	pin1: Int* = verts[(indices[prev(i, n)] & 0x0fffffff) * 4]&
	
	// If P[i] is a convex vertex [ i+1 left or on (i-1,i) ].
	if (leftOn(pin1, pi, pi1))
		return left(pi, pj, pin1) && left(pj, pi, pi1)
	// Assume (i-1,i,i+1) not collinear.
	// else P[i] is reflex.
	return !(leftOn(pi, pj, pi1) && leftOn(pj, pi, pin1))
}

// Returns T iff (v_i, v_j) is a proper internal
// diagonal of P.
diagonal: static func(i, j, n: Int, verts, indices: Int*) -> Bool {
	return inCone(i, j, n, verts, indices) && diagonalie(i, j, n, verts, indices)
}

triangulate: func(n: Int, verts, indices, tris: Int*) -> Int{
	ntris := 0
	dst := tris
	
	// The last bit of the index is used to indicate if the vertex can be removed.
	for (i: Int in 0..n) {
		i1 := next(i, n)
		i2 := next(i1, n)
		if (diagonal(i, i2, n, verts, indices))
			indices[i1] |= 0x80000000
	}
	
	while (n > 3) {
		minLen := -1
		mini := -1
		for (i: Int in 0..n) {
			i1 := next(i, n)
			if (indices[i1] & 0x80000000) {
				p0: Int* = verts[(indices[i] & 0x0fffffff) * 4]&
				p2: Int* = verts[(indices[next(i1, n)] & 0x0fffffff) * 4]&
				
				dx := p2[0] - p0[0]
				dy := p2[2] - p0[2]
				len := dx*dx + dy*dy
				
				if (minLen < 0 || len < minLen) {
					minLen = len
					mini = i
				}
			}
		}
		
		if (mini == -1) {
			// Should not happen.
			if (rcGetLog())
				rcGetLog() log(RCLogCategory RC_LOG_WARNING, "triangulate: Failed to triangulate polygon.")
/*			printf("mini == -1 ntris=%d n=%d\n", ntris, n)
			for (int i = 0; i < n; i++)
			{
				printf("%d ", indices[i] & 0x0fffffff)
			}
			printf("\n");*/
			return -ntris
		}
		
		i := mini
		i1 := next(i, n)
		i2 := next(i1, n)
		
		dst@ = indices[i] & 0x0fffffff; dst+=1
		dst@ = indices[i1] & 0x0fffffff; dst+=1
		dst@ = indices[i2] & 0x0fffffff; dst+=1
		ntris+=1
		
		// Removes P[i1] by copying P[i+1]...P[n-1] left one index.
		n-=1
		for (k: Int in i1..n)
			indices[k] = indices[k+1]
		
		if (i1 >= n) i1 = 0
		i = prev(i1, n)
		// Update diagonal flags.
		if (diagonal(prev(i, n), i1, n, verts, indices))
			indices[i] |= 0x80000000
		else
			indices[i] &= 0x0fffffff
		
		if (diagonal(i, next(i1, n), n, verts, indices))
			indices[i1] |= 0x80000000
		else
			indices[i1] &= 0x0fffffff
	}
	
	// Append the remaining triangle.
	dst = indices[0] & 0x0fffffff; dst+=1
	dst@ = indices[1] & 0x0fffffff; dst+=1
	dst@ = indices[2] & 0x0fffffff; dst+=1
	ntris+=1
	
	return ntris
}

countPolyVerts: static func(p: UShort*, nvp: Int) -> Int {
	for (i: Int in 0..nvp)
		if (p[i] == RC_MESH_NULL_IDX)
			return i
	return nvp
}

uleft: inline func(a, b, c: UShort*) -> Bool {
	return ((b[0] as Int) - (a[0] as Int)) * ((c[2] as Int) - (a[2] as Int)) -
		   ((c[0] as Int) - (a[0] as Int)) * ((b[2] as Int) - (a[2] as Int)) < 0
}

getPolyMergeValue: static func(pa, pb, verts: UShort*, ea, eb: Int@, nvp: Int) -> Int {
	na := countPolyVerts(pa, nvp)
	nb := countPolyVerts(pb, nvp)
	
	// If the merged polygon would be too big, do not merge.
	if (na+nb-2 > nvp)
		return -1
	
	// Check if the polygons share an edge.
	ea = -1
	eb = -1
	
	for (i: Int in 0..na) {
		va0: UShort = pa[i]
		va1: UShort = pa[(i+1) % na]
		if (va0 > va1)
			rcSwap(va0, va1)
		for (j: Int in 0..nb) {
			vb0: UShort = pb[j]
			vb1: UShort = pb[(j+1) % nb]
			if (vb0 > vb1)
				rcSwap(vb0, vb1)
			if (va0 == vb0 && va1 == vb1) {
				ea = i
				eb = j
				break
			}
		}
	}
	
	// No common edge, cannot merge.
	if (ea == -1 || eb == -1)
		return -1
	
	// Check to see if the merged polygon would be convex.
	va, vb, vc: UShort
	
	va = pa[(ea+na-1) % na]
	vb = pa[ea]
	vc = pb[(eb+2) % nb]
	if (!uleft(verts[va*3]&, verts[vb*3]&, verts[vc*3]&))
		return -1
	
	va = pb[(eb+nb-1) % nb]
	vb = pb[eb]
	vc = pa[(ea+2) % na]
	if (!uleft(verts[va*3]&, verts[vb*3]&, verts[vc*3]&))
		return -1
	
	va = pa[ea]
	vb = pa[(ea+1)%na]
	
	dx: Int = (verts[va*3+0] as Int) - (verts[vb*3+0] as Int)
	dy: Int = (verts[va*3+2] as Int) - (verts[vb*3+2] as Int)
	
	return dx*dx + dy*dy
}

mergePolys: static func(pa, pb: UShort*, ea, eb: Int, tmp: UShort*, nvp: Int) {
	na := countPolyVerts(pa, nvp)
	nb := countPolyVerts(pb, nvp)
	
	// Merge polygons.
	memset(tmp, 0xff, sizeof(UShort)*nvp)
	n := 0
	// Add pa
	for (i: Int in 0..(na-1))
		tmp[n] = pa[(ea+1+i) % na]; n+=1
	// Add pb
	for (i: Int in 0..(nb-1))
		tmp[n] = pb[(eb+1+i) % nb]; n+=1
	
	memcpy(pa, tmp, sizeof(UShort)*nvp)
}

pushFront: static func(v: Int, arr: Int*, an: Int@) {
	an+=1
	for (i := an-1; i > 0; i-=1)
		arr[i] = arr[i-1]
	arr[0] = v
}

pushBack: static func(v: Int, arr: Int*, an: Int@) {
	arr[an] = v
	an+=1
}

removeVertex: static func(mesh: RCPolyMesh@, rem: UShort, maxTris: Int) -> Bool {
	nvp: static Int = mesh nvp
	
	// Count number of polygons to remove.
	nrem := 0
	for (i: Int in 0..mesh npolys) {
		p: UShort* = mesh polys[i*nvp*2]&
		for (j: Int in 0..nvp)
			if (p[j] == rem) { nrem+=1; break; }
	}
	
	nedges := 0
	edges := Array<Int> new(nrem*nvp*4)
	if (!edges) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_WARNING, "removeVertex: Out of memory 'edges' (%d).", nrem*nvp*4)
		return false
	}

	nhole := 0
	hole := Array<Int> new(nrem*nvp)
	if (!hole) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_WARNING, "removeVertex: Out of memory 'hole' (%d).", nrem*nvp)
		return false
	}
	
	nhreg := 0
	hreg := Array<Int> new(nrem*nvp)
	if (!hreg) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_WARNING, "removeVertex: Out of memory 'hreg' (%d).", nrem*nvp)
		return false
	}

	nharea := 0
	harea := Array<Int> new(nrem*nvp)
	if (!harea) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_WARNING, "removeVertex: Out of memory 'harea' (%d).", nrem*nvp)
		return false
	}
	
	for (i: Int in 0..(mesh npolys)) {
		p: UShort* = mesh polys[i*nvp*2]&
		nv := countPolyVerts(p, nvp)
		hasRem := false
		for (j: Int in 0..nv)
			if (p[j] == rem) hasRem = true
		if (hasRem) {
			// Collect edges which does not touch the removed vertex.
			k := nv-1
			for (j: Int in 0..nv) {
				if (p[j] != rem && p[k] != rem) {
					e: Int* = edges[nedges*4]&
					e[0] = p[k]
					e[1] = p[j]
					e[2] = mesh regs[i]
					e[3] = mesh areas[i]
					nedges+=1
				}
				k = j
			}
			// Remove the polygon.
			p2: UShort* = mesh polys[(mesh npolys-1)*nvp*2]&
			memcpy(p, p2, sizeof(UShort)*nvp)
			mesh regs[i] = mesh regs[mesh npolys-1]
			mesh areas[i] = mesh areas[mesh npolys-1]
			mesh npolys-=1
			i-=1
		}
	}
	
	// Remove vertex.
	for (i: Int in (rem as Int)..(mesh nverts)) {
		mesh verts[i*3+0] = mesh verts[(i+1)*3+0]
		mesh verts[i*3+1] = mesh verts[(i+1)*3+1]
		mesh verts[i*3+2] = mesh verts[(i+1)*3+2]
	}
	nverts-=1
	
	// Adjust indices to match the removed vertex layout.
	for (i: Int in 0..(mesh npolys)) {
		p: UShort* = mesh polys[i*nvp*2]&
		int nv = countPolyVerts(p, nvp)
		for (j: Int in 0..nv)
			if (p[j] > rem) p[j]-=1
	}
	for (i: Int in 0..nedges) {
		if (edges[i*4+0] > rem) edges[i*4+0] -= 1
		if (edges[i*4+1] > rem) edges[i*4+1] -= 1
	}

	if (nedges == 0)
		return true

	// Start with one vertex, keep appending connected
	// segments to the start and end of the hole.
	hole[nhole] = edges[0]
	hreg[nhole] = edges[2]
	harea[nhole] = edges[3]
	nhole+=1
	
	while (nedges) {
		smatch := false
		for (i: Int in 0..nedges) {
			ea: Int = edges[i*4+0]
			eb: Int = edges[i*4+1]
			r: Int = edges[i*4+2]
			a: Int = edges[i*4+3]
			add := false
			if (hole[0] == eb) {
				// The segment matches the beginning of the hole boundary.
				pushFront(ea, hole, nhole)
				pushFront(r, hreg, nhreg)
				pushFront(a, harea, nharea)
				add = true
			} else if (hole[nhole-1] == ea) {
				// The segment matches the end of the hole boundary.
				pushBack(eb, hole, nhole)
				pushBack(r, hreg, nhreg)
				pushBack(a, harea, nharea)
				add = true
			}
			if (add) {
				// The edge segment was added, remove it.
				edges[i*4+0] = edges[(nedges-1)*4+0]
				edges[i*4+1] = edges[(nedges-1)*4+1]
				edges[i*4+2] = edges[(nedges-1)*4+2]
				edges[i*4+3] = edges[(nedges-1)*4+3]
				nedges-=1
				smatch = true
				i-=1
			}
		}
		if (!smatch)
			break
	}
	
	tris := rcAllocArray(Int, nhole*3)
	if (!tris) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_WARNING, "removeVertex: Out of memory 'tris' (%d).", nhole*3)
		return false
	}

	tverts := rcAllocArray(Int, nhole*4)
	if (!tverts) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_WARNING, "removeVertex: Out of memory 'tverts' (%d).", nhole*4)
		return false
	}

	thole := rcAllocArray(Int, nhole)
	if (!tverts) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_WARNING, "removeVertex: Out of memory 'thole' (%d).", nhole)
		return false
	}
	
	// Generate temp vertex array for triangulation.
	for (i: Int in 0..nhole) {
		pi: Int = hole[i]
		tverts[i*4+0] = mesh verts[pi*3+0]
		tverts[i*4+1] = mesh verts[pi*3+1]
		tverts[i*4+2] = mesh verts[pi*3+2]
		tverts[i*4+3] = 0
		thole[i] = i
	}

	// Triangulate the hole.
	ntris := triangulate(nhole, tverts[0]&, thole[0]&, tris)
	if (ntris < 0) {
		ntris = -ntris
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_WARNING, "removeVertex: triangulate() returned bad results.")
	}
	
	// Merge the hole triangles back to polygons.
	polys := rcAllocArray(UShort, (ntris+1)*nvp)
	if (!polys) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "removeVertex: Out of memory 'polys' (%d).", (ntris+1)*nvp)
		return false
	}
	pregs := rcAlloArray(UShort, ntris)
	if (!pregs) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "removeVertex: Out of memory 'pregs' (%d).", ntris)
		return false
	}
	pareas := rcAllocArray(UInt8, ntris)
	if (!pregs) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "removeVertex: Out of memory 'pareas' (%d).", ntris)
		return false
	}
	
	tmpPoly: UShort* = polys[ntris*nvp]&
	
	// Build initial polygons.
	npolys: Int = 0
	memset(polys, 0xff, ntris*nvp*sizeof(UShort))
	for (j: Int in 0..ntris) {
		t: Int* = tris[j*3]&
		if (t[0] != t[1] && t[0] != t[2] && t[1] != t[2]) {
			polys[npolys*nvp+0] = hole[t[0]] as UShort
			polys[npolys*nvp+1] = hole[t[1]] as UShort
			polys[npolys*nvp+2] = hole[t[2]] as UShort
			pregs[npolys] = hreg[t[0]] as UShort
			pareas[npolys] = harea[t[0]] as UInt8
			npolys+=1
		}
	}
	if (!npolys)
		return true
	
	// Merge polygons.
	if (nvp > 3) {
		while (true) {
			// Find best polygons to merge.
			bestMergeVal: Int = 0
			bestPa, bestPb, bestEa, bestEb: Int
			
			for (j: Int in 0..(npolys-1)) {
				pj: UShort* = polys[j*nvp]&
				for (k: Int in (j+1)..npolys) {
					pk: UShort* = polys[k*nvp]&
					ea, eb: Int
					v := getPolyMergeValue(pj, pk, mesh verts, ea, eb, nvp)
					if (v > bestMergeVal) {
						bestMergeVal = v
						bestPa = j
						bestPb = k
						bestEa = ea
						bestEb = eb
					}
				}
			}
			
			if (bestMergeVal > 0) {
				// Found best, merge.
				pa: UShort* = polys[bestPa*nvp]&
				pb: UShort* = polys[bestPb*nvp]&
				mergePolys(pa, pb, bestEa, bestEb, tmpPoly, nvp)
				memcpy(pb, polys[(npolys-1)*nvp]&, sizeof(UShort)*nvp)
				pregs[bestPb] = pregs[npolys-1]
				pareas[bestPb] = pareas[npolys-1]
				npolys-=1
			} else {
				// Could not merge any polygons, stop.
				break
			}
		}
	}
	
	// Store polygons.
	for (i: Int in 0..npolys) {
		if (mesh npolys >= maxTris) break
		p: UShort* = mesh polys[(mesh npolys)*nvp*2]&
		memset(p, 0xff, sizeof(UShort)*nvp*2)
		for (j: Int in 0..nvp)
			p[j] = polys[i*nvp+j]
		mesh regs[mesh npolys] = pregs[i]
		mesh areas[mesh npolys] = pareas[i]
		mesh npolys+=1
		if (mesh npolys > maxTris) {
			if (rcGetLog())
				rcGetLog() log(RC_LOG_ERROR, "removeVertex: Too many polygons %d (max:%d).", mesh npolys, maxTris)
			return false
		}
	}
	return true
}

rcBuildPolyMesh: func(cset: RCContourSet@, nvp: Int, mesh: RCPolyMesh) -> Bool {
	startTime := rcGetPerformanceTimer()
	
	vcopy(mesh bmin, cset bmin)
	vcopy(mesh bmax, cset bmax)
	mesh cs = cset cs
	mesh ch = cset ch
	
	maxVertices := 0
	maxTris := 0
	maxVertsPerCont := 0
	for (i: Int in 0..(cset nconts)) {
		// Skip null contours.
		if (cset conts[i] nverts < 3) continue
		maxVertices += cset conts[i] nverts
		maxTris += cset conts[i] nverts - 2
		maxVertsPerCont = rcMax(maxVertsPerCont, cset conts[i] nverts)
	}
	
	if (maxVertices >= 0xfffe) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMesh: Too many vertices %d.", maxVertices)
		return false
	}
	
	vflags := Array<UInt8> new(maxVertices)
	if (!vflags) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMesh: Out of memory 'mesh verts' (%d).", maxVertices)
		return false
	}
	memset(vflags, 0, maxVertices)
	
	mesh verts = Array<UShort> new(maxVertices*3)
	if (!mesh verts) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMesh: Out of memory 'mesh verts' (%d).", maxVertices)
		return false
	}
	
	mesh polys = Array<UShort> new(maxTris*nvp*2*2)
	if (!mesh polys) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMesh: Out of memory 'mesh polys' (%d).", maxTris*nvp*2)
		return false
	}
	
	mesh regs = Array<UShort> new(maxTris)
	if (!mesh regs) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMesh: Out of memory 'mesh.regs' (%d).", maxTris)
		return false
	}
	
	mesh areas = Array<UInt8> new(maxTris)
	if (!mesh areas) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMesh: Out of memory 'mesh.areas' (%d).", maxTris)
		return false
	}
	
	mesh nverts = 0
	mesh npolys = 0
	mesh nvp = nvp
	
	memset(mesh verts, 0, sizeof(UShort)*maxVertices*3)
	memset(mesh polys, 0xff, sizeof(UShort)*maxTris*nvp*2)
	memset(mesh regs, 0, sizeof(UShort)*maxTris)
	memset(mesh areas, 0, sizeof(UInt8)*maxTris)
	
	nextVert := Array<Int> new(maxVertices)
	if (!nextVert) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMesh: Out of memory 'nextVert' (%d).", maxVertices)
		return false
	}
	memset(nextVert, 0, sizeof(Int)*maxVertices)
	
	firstVert := Array<Int> new(VERTEX_BUCKET_COUNT)
	if (!firstVert) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMesh: Out of memory 'firstVert' (%d).", VERTEX_BUCKET_COUNT)
		return false
	}
	for (i: Int in 0..VERTEX_BUCKET_COUNT)
		firstVert[i] = -1
	
	indices := Array<Int> new(maxVertsPerCont)
	if (!indices) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMesh: Out of memory 'indices' (%d).", maxVertsPerCont)
		return false
	}
	
	tris := Array<Int> new(maxVertsPerCont*3)
	if (!tris) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMesh: Out of memory 'tris' (%d).", maxVertsPerCont*3)
		return false
	}
	
	polys := Array<UShort> new((maxVertsPerCont+1)*nvp) 
	if (!polys) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMesh: Out of memory 'polys' (%d).", (maxVertsPerCont+1)*nvp)
		return false
	}
	tmpPoly: UShort* = polys[maxVertsPerCont*nvp]&
	
	for (i: Int in 0..(cset nconts)) {
		cont: RCContour = cset conts[i]&
		
		// Skip null contours.
		if (cont nverts < 3)
			continue
		
		// Triangulate contour
		for (j: Int in 0..(cont nverts))
			indices[j] = j
			
		ntris := triangulate(cont nverts, cont verts, indices[0]&, tris[0]&)
		if (ntris <= 0) {
			// Bad triangulation, should not happen.
			/*for (int k = 0; k < cont.nverts; ++k)
			{
				int* v = cont verts[k*4]&
				printf("\t\t%d,%d,%d,%d,\n", v[0], v[1], v[2], v[3])
				if (nBadPos < 100)
				{
					badPos[nBadPos*3+0] = v[0]
					badPos[nBadPos*3+1] = v[1]
					badPos[nBadPos*3+2] = v[2]
					nBadPos++
				}
			}*/
			ntris = -ntris
		}
		
		// Add and merge vertices.
		for (j: Int in 0..(cont nverts)) {
			v: Int* = cont verts[j*4]&
			indices[j] = addVertex(v[0] as UShort, v[1] as UShort, v[2] as UShort,
								   mesh verts, firstVert, nextVert, mesh nverts)
			if (v[3] & RC_BORDER_VERTEX) {
				// This vertex should be removed.
				vflags[indices[j]] = 1
			}
		}
		
		// Build initial polygons.
		npolys := 0
		memset(polys, 0xff, maxVertsPerCont*nvp*sizeof(UShort))
		for (j: Int in 0..ntris) {
			int* t = tris[j*3]&
			if (t[0] != t[1] && t[0] != t[2] && t[1] != t[2]) {
				polys[npolys*nvp+0] = indices[t[0]] as UShort
				polys[npolys*nvp+1] = indices[t[1]] as UShort
				polys[npolys*nvp+2] = indices[t[2]] as UShort
				npolys+=1
			}
		}
		if (!npolys)
			continue
		
		// Merge polygons.
		if (nvp > 3) {
			while (true) {
				// Find best polygons to merge.
				bestMergeVal := 0
				bestPa, bestPb, bestEa, bestEb: Int
				
				for (j: Int in 0..npolys-1) {
					pj: UShort* = polys[j*nvp]&
					for (k in j+1..npolys) {
						UShort* pk = polys[k*nvp]&
						ea, eb: Int
						v := getPolyMergeValue(pj, pk, mesh verts, ea, eb, nvp)
						if (v > bestMergeVal) {
							bestMergeVal = v
							bestPa = j
							bestPb = k
							bestEa = ea
							bestEb = eb
						}
					}
				}
				
				if (bestMergeVal > 0) {
					// Found best, merge.
					UShort* pa = polys[bestPa*nvp]&
					UShort* pb = polys[bestPb*nvp]&
					mergePolys(pa, pb, bestEa, bestEb, tmpPoly, nvp)
					memcpy(pb, polys[(npolys-1)*nvp]&, sizeof(UShort)*nvp)
					npolys-=1
				} else {
					// Could not merge any polygons, stop.
					break
				}
			}
		}
		
		// Store polygons.
		for (j: Int in 0..npolys) {
			p: UShort* = mesh polys[mesh npolys*nvp*2]&
			q: UShort* = polys[j*nvp]&
			for (k in 0..nvp)
				p[k] = q[k]
			mesh regs[mesh npolys] = cont reg
			mesh areas[mesh npolys] = cont area
			mesh npolys+=1
			if (mesh npolys > maxTris) {
				if (rcGetLog())
					rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMesh: Too many polygons %d (max:%d).", mesh npolys, maxTris)
				return false
			}
		}
	}
	
	// Remove edge vertices.
	for (i: Int in 0..(mesh nverts)) {
		if (vflags[i]) {
			if (!removeVertex(mesh, i, maxTris)) {
				if (rcGetLog())
					rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMesh: Failed to remove edge vertex %d.", i)
				return false
			}
			// Note: mesh nverts is already decremented inside removeVertex()!
			for (j: Int in i..(mesh nverts))
				vflags[j] = vflags[j+1]
			--i
		}
	}
	
	// Calculate adjacency.
	if (!buildMeshAdjacency(mesh polys, mesh npolys, mesh nverts, nvp)) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMesh: Adjacency failed.")
		return false
	}
	
	endTime := rcGetPerformanceTimer()
	
//	if (rcGetLog())
//		rcGetLog() log(RC_LOG_PROGRESS, "Build polymesh: %.3f ms", rcGetDeltaTimeUsec(startTime, endTime)/1000.0f)
	if (rcGetBuildTimes())
		rcGetBuildTimes() buildPolymesh += rcGetDeltaTimeUsec(startTime, endTime)
	
	return true
}

rcMergePolyMeshes: func(meshes: RCPolyMesh**, nmeshes: Int, mesh: RCPolyMesh@) -> Bool {
	if (!nmeshes || !meshes)
		return true
	
	startTime := rcGetPerformanceTimer()
	
	mesh nvp = meshes[0] nvp
	mesh cs = meshes[0] cs
	mesh ch = meshes[0] ch
	vcopy(mesh bmin, meshes[0] bmin)
	vcopy(mesh bmax, meshes[0] bmax)
	
	maxVerts := 0
	maxPolys := 0
	maxVertsPerMesh := 0
	for (i: Int in 0..nmeshes) {
		vmin(mesh bmin, meshes[i] bmin)
		vmax(mesh bmax, meshes[i] bmax)
		maxVertsPerMesh = rcMax(maxVertsPerMesh, meshes[i] nverts)
		maxVerts += meshes[i] nverts
		maxPolys += meshes[i] npolys
	}
	
	mesh nverts = 0
	mesh verts = rcAllocArray(UShort, maxVerts*3)
	if (!mesh verts) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcMergePolyMeshes: Out of memory 'mesh verts' (%d).", maxVerts*3)
		return false
	}
	
	mesh npolys = 0
	mesh polys = rcAllocArray(UShort, maxPolys*2*mesh nvp)
	if (!mesh polys) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcMergePolyMeshes: Out of memory 'mesh polys' (%d).", maxPolys*2*mesh nvp)
		return false
	}
	memset(mesh polys, 0xff, sizeof(UShort)*maxPolys*2*mesh nvp)
	
	mesh regs = rcAllocArray(UShort, maxPolys)
	if (!mesh regs) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcMergePolyMeshes: Out of memory 'mesh.regs' (%d).", maxPolys)
		return false
	}
	memset(mesh regs, 0, sizeof(UShort)*maxPolys)
	
	mesh areas = rcAllocArray(UInt8, maxPolys)
	if (!mesh areas) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcMergePolyMeshes: Out of memory 'mesh.areas' (%d).", maxPolys)
		return false
	}
	memset(mesh areas, 0, sizeof(UInt8)*maxPolys)
	
	nextVert := rcAllocArray(Int, maxVerts)
	if (!nextVert) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcMergePolyMeshes: Out of memory 'nextVert' (%d).", maxVerts)
		return false
	}
	memset(nextVert, 0, sizeof(Int)*maxVerts)
	
	firstVert := rcAllocArray(Int, VERTEX_BUCKET_COUNT)
	if (!firstVert) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcMergePolyMeshes: Out of memory 'firstVert' (%d).", VERTEX_BUCKET_COUNT)
		return false
	}
	
	for (i: Int in 0..VERTEX_BUCKET_COUNT)
		firstVert[i] = -1
	
	vremap := rcAllocArray(UShort, maxVertsPerMesh)
	if (!vremap) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcMergePolyMeshes: Out of memory 'vremap' (%d).", maxVertsPerMesh)
		return false
	}
	memset(nextVert, 0, sizeof(Int)*maxVerts)
	
	for (i: Int in 0..nmeshes) {
		pmesh: RCPolyMesh* = meshes[i]
		
		ox: UShort = floor((pmesh bmin[0]-mesh bmin[0])/mesh cs+0.5) as UShort
		oz: UShort = floor((pmesh bmin[2]-mesh bmin[2])/mesh cs+0.5) as UShort
		
		for (j: Int in 0..(pmesh nverts)) {
			v: UShort* = pmesh verts[j*3]&
			vremap[j] = addVertex(v[0]+ox, v[1], v[2]+oz,
								  mesh verts, firstVert, nextVert, mesh nverts)
		}
		
		for (j: Int in 0..(pmesh npolys)) {
			tgt: UShort* = mesh polys[mesh npolys*2*mesh nvp]&
			src: UShort* = pmesh polys[j*2*mesh nvp]&
			mesh regs[mesh npolys] = pmesh regs[j]
			mesh areas[mesh npolys] = pmesh areas[j]
			mesh npolys+=1
			for (k in 0..(mesh nvp)) {
				if (src[k] == RC_MESH_NULL_IDX) break
				tgt[k] = vremap[src[k]]
			}
		}
	}
	
	// Calculate adjacency.
	if (!buildMeshAdjacency(mesh polys, mesh npolys, mesh nverts, mesh nvp)) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcMergePolyMeshes: Adjacency failed.")
		return false
	}
	
	endTime := rcGetPerformanceTimer()
	if (rcGetBuildTimes())
		rcGetBuildTimes() mergePolyMesh += rcGetDeltaTimeUsec(startTime, endTime)
	
	return true
}

