use retour
import Retour/[Retour]

RC_UNSET_HEIGHT: static const UShort = 0xffff

RCHeightPatch: class {
	init: func() {}
	destroy: func() { delete [] data; }
	data: UShort*
	xmin, ymin, width, height: Int
}

vdot2: inline func(a, b: Float*) -> Float {
	return a[0]*b[0] + a[2]*b[2]
}

vdistSq2: inline func(p, q: Float*) -> Float {
	dx := q[0] - p[0]
	dy := q[2] - p[2]
	return dx*dx + dy*dy
}

vdist2: inline func(p, q: Float*) -> Float {
	return sqrt(vdistSq2(p, q))
}

vcross2: inline func(p1, p2, p3: Float*) -> Float {
	u1 := p2[0] - p1[0]
	v1 := p2[2] - p1[2]
	u2 := p3[0] - p1[0]
	v2 := p3[2] - p1[2]
	return u1 * v2 - v1 * u2
}

circumCircle: static func(p1, p2, p3, c: Float*, r: Float@) -> Bool {
	EPS: static const Float = 1e-6
	
	cp := vcross2(p1, p2, p3)
	if (cp abs() > EPS) {
		p1Sq := vdot2(p1, p1)
		p2Sq := vdot2(p2, p2)
		p3Sq := vdot2(p3, p3)
		c[0] = (p1Sq*(p2[2]-p3[2]) + p2Sq*(p3[2]-p1[2]) + p3Sq*(p1[2]-p2[2])) / (2*cp)
		c[2] = (p1Sq*(p3[0]-p2[0]) + p2Sq*(p1[0]-p3[0]) + p3Sq*(p2[0]-p1[0])) / (2*cp)
		r = vdist2(c, p1)
		return true
	}
	
	c[0] = p1[0]
	c[2] = p1[2]
	r = 0
	return false
}

distPtTri: static func(p, a, b, c: Float*) -> Float {
	EPS: static const Float = 1e-4
	
	v0, v1, v2: Float[3]
	vsub(v0, c, a)
	vsub(v1, b, a)
	vsub(v2, p, a)
	
	dot00 := vdot2(v0, v0)
	dot01 := vdot2(v0, v1)
	dot02 := vdot2(v0, v2)
	dot11 := vdot2(v1, v1)
	dot12 := vdot2(v1, v2)
	
	// Compute barycentric coordinates
	invDenom: Float = 1.0 / (dot00 * dot11 - dot01 * dot01)
	u: Float = (dot11 * dot02 - dot01 * dot12) * invDenom
	v: Float = (dot00 * dot12 - dot01 * dot02) * invDenom
	
	// If poInt lies inside the triangle, return Interpolated y-coord.
	if (u >= -EPS && v >= -EPS && (u+v) <= 1+EPS) {
		y := a[1] + v0[1]*u + v1[1]*v
		return (y-p[1]) abs()
	}
	return FLT_MAX
}

distancePtSeg: static func(pt, p, q: Float*) -> Float {
	pqx := q[0] - p[0]
	pqy := q[1] - p[1]
	pqz := q[2] - p[2]
	dx := pt[0] - p[0]
	dy := pt[1] - p[1]
	dz := pt[2] - p[2]
	d := pqx*pqx + pqy*pqy + pqz*pqz
	t := pqx*dx + pqy*dy + pqz*dz
	if (d > 0)
		t /= d
	if (t < 0)
		t = 0
	else if (t > 1)
		t = 1
	
	dx = p[0] + t*pqx - pt[0]
	dy = p[1] + t*pqy - pt[1]
	dz = p[2] + t*pqz - pt[2]
	
	return dx*dx + dy*dy + dz*dz
}

distancePtSeg2d: static func(pt, p, q: Float*) -> Float {
	pqx := q[0] - p[0]
	pqz := q[2] - p[2]
	dx := pt[0] - p[0]
	dz := pt[2] - p[2]
	d := pqx*pqx + pqz*pqz
	t := pqx*dx + pqz*dz
	if (d > 0)
		t /= d
	if (t < 0)
		t = 0
	else if (t > 1)
		t = 1
	
	dx = p[0] + t*pqx - pt[0]
	dz = p[2] + t*pqz - pt[2]
	
	return dx*dx + dz*dz
}

distToTriMesh: static func(p, verts: Float*, nverts: Int, tris: Int*, ntris: Int) -> Float {
	dmin := FLT_MAX
	for (i: Int in 0. ntris) {
		va: Float* = verts[tris[i*4+0]*3]&
		vb: Float* = verts[tris[i*4+1]*3]&
		vc: Float* = verts[tris[i*4+2]*3]&
		d := distPtTri(p, va, vb, vc)
		if (d < dmin)
			dmin = d
	}
	if (dmin == FLT_MAX) return -1
	return dmin
}

distToPoly: static func(nvert: Int, verts, p: Float*) -> Float {
	dmin := FLT_MAX
	c, j: Int
	j = nvert-1
	for (i: Int in 0..nvert) {
		vi: Float* = verts[i*3]&
		vj: Float* = verts[j*3]&
		if (((vi[2] > p[2]) != (vj[2] > p[2])) &&
			(p[0] < (vj[0]-vi[0]) * (p[2]-vi[2]) / (vj[2]-vi[2]) + vi[0]))
			c = !c
		dmin = rcMin(dmin, distancePtSeg2d(p, vj, vi))
		j = i
	}
	return c ? -dmin : dmin
}

getHeight: static func(fx, fz, cs, ics: Float, hp: RCHeightPatch@) -> UShort {
	ix := floor(fx*ics + 0.01) as Int
	iz := floor(fz*ics + 0.01) as Int
	ix = rcClamp(ix-hp xmin, 0, hp width)
	iz = rcClamp(iz-hp ymin, 0, hp height)
	h: UShort = hp data[ix+iz*hp width]
	if (h == RC_UNSET_HEIGHT) {
		// Special case when data might be bad.
		// Find nearest neighbour pixel which has valid height.
		off: Int[16] = [-1,0, -1,-1, 0,-1, 1,-1, 1,0, 1,1, 0,1, -1,1]
		dmin := FLT_MAX
		for (i: Int in 0..8) {
			nx := ix+off[i*2+0]
			nz := iz+off[i*2+1]
			if (nx < 0 || nz < 0 || nx >= hp width || nz >= hp height) continue
			nh := hp data[nx+nz*hp width]
			if (nh == RC_UNSET_HEIGHT) continue
			dx := (nx+0.5)*cs - fx
			dz := (nz+0.5)*cs - fz
			d := dx*dx+dz*dz
			if (d < dmin) {
				h = nh
				dmin = d
			} 
		}
	}
	return h
}

//enum
EdgeValues: class {
	UNDEF: static const Int = -1
	HULL: static const Int = -2
}

findEdge: static func(edges: Int*, nedges, s, t: Int) -> Int {
	for (i: Int in 0..nedges) {
		e: Int* = edges[i*4]&
		if ((e[0] == s && e[1] == t) || (e[0] == t && e[1] == s))
			return i
	}
	return EdgeValues UNDEF
}

addEdge: static func(edges: Int*, nedges: Int@, maxEdges, s, t, l, r: Int) -> Int {
	if (nedges >= maxEdges) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "addEdge: Too many edges (%d/%d).", nedges, maxEdges)
		return EdgeValues UNDEF
	}
	
	// Add edge if not already in the triangulation. 
	e := findEdge(edges, nedges, s, t)
	if (e == EdgeValues UNDEF) {
		e: Int* = edges[nedges*4]&
		e[0] = s
		e[1] = t
		e[2] = l
		e[3] = r
		sretval := nedges
		nedges+=1
		return sretval
	} else {
		return EdgeValues UNDEF
	}
}

updateLeftFace: static func(e: Int*, s, t, f: Int) {
	if (e[0] == s && e[1] == t && e[2] == EdgeValues UNDEF)
		e[2] = f
	else if (e[1] == s && e[0] == t && e[3] == EdgeValues UNDEF)
		e[3] = f
}	

overlapSegSeg2d: static func(a, b, c, d: Float*) -> Int {
	a1 := vcross2(a, b, d)
	a2 := vcross2(a, b, c)
	if (a1 * a2 < 0.0) {
		a3 := vcross2(c, d, a)
		a4 := a3 + a2 - a1
		if (a3 * a4 < 0.0)
			return 1
	}	
	return 0
}

overlapEdges: static func(pts: Float*, edges: Int*, nedges, s1, t1: Int) -> Bool {
	for (i: Int in 0..nedges) {
		s0 := edges[i*4+0]
		t0 := edges[i*4+1]
		// Same or connected edges do not overlap.
		if (s0 == s1 || s0 == t1 || t0 == s1 || t0 == t1)
			continue
		if (overlapSegSeg2d(pts[s0*3]&, pts[t0*3]&, pts[s1*3]&, pts[t1*3]&))
			return true
	}
	return false
}

completeFacet: static func(pts: Float*, npts: Int, edges: Int*, nedges: Int@, maxEdges: Int, nfaces: Int@, e: Int) {
	EPS: static const Float = 1e-5
	
	edge: Int* = edges[e*4]&
	
	// Cache s and t.
	s, t: Int
	if (edge[2] == EdgeValues UNDEF) {
		s = edge[0]
		t = edge[1]
	} else if (edge[3] == EdgeValues UNDEF) {
		s = edge[1]
		t = edge[0]
	} else {
	    // Edge already completed. 
	    return
	}
    
	// Find best poInt on left of edge. 
	pt := npts
	c: Float[3] // = [0, 0, 0]
	r: Float = -1.0
	for (u: Int in 0..npts) {
		if (u == s || u == t) continue
		if (vcross2(pts[s*3]@, pts[t*3]@, pts[u*3]&) > EPS) {
			if (r < 0) {
				// The circle is not updated yet, do it now.
				pt = u
				circumCircle(pts[s*3]&, pts[t*3]&, pts[u*3]&, c, r)
				continue
			}
			d := vdist2(c, pts[u*3])
			tol: Float = 0.001
			if (d > r*(1+tol)) {
				// Outside current circumcircle, skip.
				continue
			} else if (d < r*(1-tol)) {
				// Inside safe circumcircle, update circle.
				pt = u
				circumCircle(pts[s*3]&, pts[t*3]&, pts[u*3]&, c, r)
			} else {
				// Inside epsilon circum circle, do extra tests to make sure the edge is valid.
				// s-u and t-u cannot overlap with s-pt nor t-pt if they exists.
				if (overlapEdges(pts, edges, nedges, s, u))
					continue
				if (overlapEdges(pts, edges, nedges, t, u))
					continue
				// Edge is valid.
				pt = u
				circumCircle(pts[s*3]&, pts[t*3]&, pts[u*3]&, c, r)
			}
		}
	}
	
	// Add new triangle or update edge info if s-t is on hull. 
	if (pt < npts) {
		// Update face information of edge being completed. 
		updateLeftFace(edges[e*4]&, s, t, nfaces)
		
		// Add new edge or update face info of old edge. 
		e = findEdge(edges, nedges, pt, s)
		if (e == EdgeValues UNDEF)
		    addEdge(edges, nedges, maxEdges, pt, s, nfaces, EdgeValues UNDEF)
		else
		    updateLeftFace(edges[e*4]&, pt, s, nfaces)
		
		// Add new edge or update face info of old edge. 
		e = findEdge(edges, nedges, t, pt)
		if (e == EdgeValues UNDEF)
		    addEdge(edges, nedges, maxEdges, t, pt, nfaces, EdgeValues UNDEF)
		else
		    updateLeftFace(edges[e*4]&, t, pt, nfaces)
		
		nfaces+=1
	} else {
		updateLeftFace(edges[e*4]&, s, t, EdgeValues HULL)
	}
}

delaunayHull: static func(npts: Int, pts: Float*, nhull: Int, hull: Int*, tris: RCIntArray@, edges: RCIntArray@) {
	nfaces := 0
	nedges := 0
	maxEdges := npts*10
	edges resize(maxEdges*4)
	
	j := nhull-1
	for (i: Int in 0..nhull) {
		addEdge(edges[0]&, nedges, maxEdges, hull[j], hull[i], EdgeValues HULL, EdgeValues UNDEF)
		j = i
	}
	
	currentEdge := 0
	while (currentEdge < nedges) {
		if (edges[currentEdge*4+2] == EdgeValues UNDEF)
			completeFacet(pts, npts, edges[0]&, nedges, maxEdges, nfaces, currentEdge)
		if (edges[currentEdge*4+3] == EdgeValues UNDEF)
			completeFacet(pts, npts, edges[0]&, nedges, maxEdges, nfaces, currentEdge)
		currentEdge+=1
	}

	// Create tris
	tris resize(nfaces*4)
	for (i: Int in 0..(nfaces*4))
		tris[i] = -1
	
	for (i: Int in 0..nedges) {
		e: Int* = edges[i*4]&
		if (e[3] >= 0) {
			// Left face
			t: Int* = tris[e[3]*4]&
			if (t[0] == -1) {
				t[0] = e[0]
				t[1] = e[1]
			} else if (t[0] == e[1])
				t[2] = e[0]
			else if (t[1] == e[0])
				t[2] = e[1]
		}
		if (e[2] >= 0) {
			// Right
			t: Int* = tris[e[2]*4]&
			if (t[0] == -1) {
				t[0] = e[1]
				t[1] = e[0]
			} else if (t[0] == e[0])
				t[2] = e[1]
			else if (t[1] == e[1])
				t[2] = e[0]
		}
	}
	
	for (i: Int in 0..(tris size()/4)) {
		t: Int* = tris[i*4]&
		if (t[0] == -1 || t[1] == -1 || t[2] == -1) {
			if (rcGetLog())
				rcGetLog() log(RCLogCategory RC_LOG_WARNING, "delaunayHull: Removing dangling face %d [%d,%d,%d].", i, t[0], t[1], t[2])
			t[0] = tris[tris size()-4]
			t[1] = tris[tris size()-3]
			t[2] = tris[tris size()-2]
			t[3] = tris[tris size()-1]
			tris resize(tris size()-4)
		}
	}

}



buildPolyDetail: static func(zin: Float*, nin: Int, sampleDist, sampleMaxError: Float,
							chf: RCCompactHeightfield@, hp: RCHeightPatch@, verts: Float*,
							nverts: Int@, tris, edges, samples: RCIntArray@) -> Bool {
	MAX_VERTS: static const Int = 256
	MAX_EDGE: static const Int = 64
	edge: Float[(MAX_EDGE+1)*3]
	hull: Int[MAX_VERTS]
	nhull = 0
	
	nverts = 0
	for (i: Int in 0..nin)
		vcopy(verts[i*3]&, zin[i*3]&)
	nverts = nin
	
	cs := chf cs
	ics: Float = 1.0/cs
	
	// Tesselate outlines.
	// This is done in separate pass in order to ensure
	// seamless height values across the ply boundaries.
	if (sampleDist > 0) {
		j := nin-1
		for (i: Int in 0..nin) {
			vj: Float* = zin[j*3]&
			vi: Float* = zin[i*3]&
			swapped := false
			// Make sure the segments are always handled in same order
			// using lexological sort or else there will be seams.
			if (fabsf(vj[0]-vi[0]) < 1e-6) {
				if (vj[2] > vi[2]) {
					rcSwap(vj,vi)
					swapped = true
				}
			} else {
				if (vj[0] > vi[0]) {
					rcSwap(vj,vi)
					swapped = true
				}
			}
			// Create samples along the edge.
			dx := vi[0] - vj[0]
			dz := vi[2] - vj[2]
			d := sqrtf(dx*dx + dz*dz)
			nn := 1 + floor(d/sampleDist) as Int
			if (nn > MAX_EDGE) nn = MAX_EDGE
			if (nverts+nn >= MAX_VERTS)
				nn = MAX_VERTS-1-nverts
			for (k := 0; k <= nn; k+=1) {
				u: Float = (k as Float)/(nn as Float)
				pos: Float* = edge[k*3]&
				pos[0] = vj[0] + dx*u
				pos[2] = vj[2] + dz*u
				pos[1] = getHeight(pos[0], pos[2], cs, ics, hp)*chf ch
			}
			// Simplify samples.
			idx: Int[MAX_EDGE] = [0, nn]
			nidx := 2
			k := 0
			while (k < (nidx-1)) {
				a := idx[k]
				b := idx[k+1]
				va: Float* = edge[a*3]&
				vb: Float* = edge[b*3]&
				// Find maximum deviation along the segment.
				maxd: Float = 0.0
				maxi := -1
				for (m: Int in (a+1)..b) {
					d := distancePtSeg(edge[m*3]&, va, vb)
					if (d > maxd) {
						maxd = d
						maxi = m
					}
				}
				// If the max deviation is larger than accepted error,
				// add new poInt, else continue to next segment.
				if (maxi != -1 && maxd > rcSqr(sampleMaxError)) {
					for (m := nidx; m > k; --m)
						idx[m] = idx[m-1]
					idx[k+1] = maxi
					nidx+=1
				} else {
					k+=1
				}
			}
			
			hull[nhull] = j
			nhull+=1
			// Add new vertices.
			if (swapped) {
				for (k := nidx-2; k > 0; k-=1) {
					vcopy(verts[nverts*3]&, edge[idx[k]*3]&)
					hull[nhull] = nverts
					nhull+=1; nhull+=1
				}
			} else {
				for (k: Int in 1..(nidx-1)) {
					vcopy(verts[nverts*3]&, edge[idx[k]*3]&)
					hull[nhull] = nverts
					nhull+=1; nverts+=1
				}
			}
			j=i
		}
	}
	

	// Tesselate the base mesh.
	edges resize(0)
	tris resize(0)
	delaunayHull(nverts, verts, nhull, hull, tris, edges)
	
	if (tris size() == 0) {
		// Could not triangulate the poly, make sure there is some valid data there.
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_WARNING, "buildPolyDetail: Could not triangulate polygon, adding default data.")
		for (i: Int in 2. nverts) {
			tris push(0)
			tris push(i-1)
			tris push(i)
			tris push(0)
		}
		return true
	}

	if (sampleDist > 0) {
		// Create sample locations in a grid.
		bmin, bmax: Float[3]
		vcopy(bmin, zin)
		vcopy(bmax, zin)
		for (i: Int in 1..nin) {
			vmin(bmin, zin[i*3]&)
			vmax(bmax, zin[i*3]&)
		}
		x0 := floor(bmin[0]/sampleDist) as Int
		x1 := ceil(bmax[0]/sampleDist) as Int
		z0 := floor(bmin[2]/sampleDist) as Int
		z1 := ceil(bmax[2]/sampleDist) as Int
		samples.resize(0)
		for (z: Int in z0..z1) {
			for (x: Int in x0..x1) {
				pt: Float[3]
				pt[0] = x*sampleDist
				pt[2] = z*sampleDist
				// Make sure the samples are not too close to the edges.
				if (distToPoly(nin, zin, pt) > -sampleDist/2) continue
				samples push(x)
				samples push(getHeight(pt[0], pt[2], cs, ics, hp))
				samples push(z)
			}
		}
				
		// Add the samples starting from the one that has the most
		// error. The procedure stops when all samples are added
		// or when the max error is within treshold.
		nsamples := samples size()/3
		for (iten: Int in 0..nsamples) {
			// Find sample with most error.
			bestpt: Float[3]
			bestd: Float = 0.0
			for (i: Int in 0..nsamples) {
				pt: Float[3]
				pt[0] = samples[i*3+0]*sampleDist
				pt[1] = samples[i*3+1]*chf ch
				pt[2] = samples[i*3+2]*sampleDist
				d := distToTriMesh(pt, verts, nverts, tris[0]&, tris size()/4)
				if (d < 0) continue; // did not hit the mesh.
				if (d > bestd) {
					bestd = d
					vcopy(bestpt, pt)
				}
			}
			// If the max error is within accepted threshold, stop tesselating.
			if (bestd <= sampleMaxError)
				break

			// Add the new sample poInt.
			vcopy(verts[nverts*3]&, bestpt)
			nverts+=1
			
			// Create new triangulation.
			// TODO: Incremental add instead of full rebuild.
			edges resize(0)
			tris resize(0)
			delaunayHull(nverts, verts, nhull, hull, tris, edges)
			
			if (nverts >= MAX_VERTS)
				break
		}
	}
	return true
}

getHeightData: static func(chf: RCCompactHeightfield@, poly: UShort*, npoly: Int, verts: UShort*, hp: RCHeightPatch@, stack: RCIntArray@) {
	// Floodfill the heightfield to get 2D height data,
	// starting at vertex locations as seeds.
	
	memset(hp data, 0xff, sizeof(UShort)*hp width*hp height)
	stack resize(0)
	
	// Use poly vertices as seed points for the flood fill.
	for (j: Int in 0..npoly) {
		ax := verts[poly[j]*3+0] as Int
		ay := verts[poly[j]*3+1] as Int
		az := verts[poly[j]*3+2] as Int
		if (ax < hp xmin || ax >= hp xmin+hp width ||
			az < hp ymin || az >= hp ymin+hp height)
			continue
			
		c := chf cells[ax+az*chf width]
		dmin: Int = RC_UNSET_HEIGHT
		ai := -1
		for (i: Int in (c index as Int)..(c index+c count) as Int) {
			s := chf spans[i]
			d := rcAbs(ay - (s y as Int))
			if (d < dmin) {
				ai = i
				dmin = d
			}
		}
		if (ai != -1) {
			stack push(ax)
			stack push(az)
			stack push(ai)
		}
	}
	
	// Not no match, try polygon center.
	if (stack size() == 0) {
		cx, cy, cz: Int
		for (j: Int in 0..npoly) {
			cx += verts[poly[j]*3+0] as Int
			cy += verts[poly[j]*3+1] as Int
			cz += verts[poly[j]*3+2] as Int
		}
		cx /= npoly
		cy /= npoly
		cz /= npoly
		
		if (cx >= hp xmin && cx < hp xmin+hp width &&
			cz >= hp ymin && cz < hp ymin+hp height) {
			c: RCCompactCell& = chf cells[cx+cz*chf width]
			dmin: Int = RC_UNSET_HEIGHT
			ci := -1
			for (i: Int in (c index as Int)..(c index+c count) as Int) {
				s := chf spans[i]
				d: Int = rcAbs(cy - (s y as Int))
				if (d < dmin) {
					ci = i
					dmin = d
				}
			}
			if (ci != -1) {
				stack push(cx)
				stack push(cz)
				stack push(ci)
			}
		}
	}
	
	while (stack size() > 0) {
		ci := stack pop()
		cy := stack pop()
		cx := stack pop()
		
		// Skip already visited locations.
		idx := cx-hp xmin+(cy-hp ymin)*hp width
		if (hp data[idx] != RC_UNSET_HEIGHT)
			continue
		
		cs := chf spans[ci]
		hp data[idx] = cs y
		
		for (dir: Int in 0..4) {
			if (rcGetCon(cs, dir) == RC_NOT_CONNECTED) continue
			
			ax := cx + rcGetDirOffsetX(dir)
			ay := cy + rcGetDirOffsetY(dir)
			if (ax < hp xmin || ax >= (hp xmin+hp width) ||
				ay < hp ymin || ay >= (hp ymin+hp height))
				continue

			if (hp data[ax-hp xmin+(ay-hp ymin)*hp width] != RC_UNSET_HEIGHT)
				continue

			ai = (chf cells[ax+ay*chf width] index as Int) + rcGetCon(cs, dir)
			stack push(ax)
			stack push(ay)
			stack push(ai)
		}
	}	
}

getEdgeFlags: static func(va, vb, vpoly: Float*, npoly: Int) -> UInt8 {
	// Return true if edge (va, vb) is part of the polygon.
	thrSqr: static Float = rcSqr(0.001)
	j: Int = npoly-1
	for (i: Int in 0..npoly) {
		if (distancePtSeg2d(va, vpoly[j*3]&, vpoly[i*3]&) < thrSqr && 
			distancePtSeg2d(vb, vpoly[j*3]&, vpoly[i*3]&) < thrSqr)
			return 1
		j = i
	}
	return 0
}

getTriFlags: static func(va, vb, vc, vpoly: Float*, npoly: Int) -> UInt8 {
	flags: UInt8 = 0
	flags |= getEdgeFlags(va, vb, vpoly, npoly) << 0
	flags |= getEdgeFlags(vb, vc, vpoly, npoly) << 2
	flags |= getEdgeFlags(vc, va, vpoly, npoly) << 4
	return flags
}

rcBuildPolyMeshDetail: func(mesh: RCPolyMesh@, chf: RCCompactHeightfield@, sampleDist, sampleMaxError: Float, dmesh: RCPolyMeshDetail@) -> Bool {
	startTime := rcGetPerformanceTimer()
	
	if (mesh nverts == 0 || mesh npolys == 0)
		return true
	
	nvp := mesh nvp
	cs := mesh cs
	ch := mesh ch
	orig: Float* = mesh bmin
	
	edges := RCIntArray new(64)
	tris := RCIntArray new(512)
	stack := RCIntArray new(512)
	samples := RCIntArray new(512)
	verts: Float[256*3]
	hp: RCHeightPatch
	nPolyVerts := 0
	maxhw := 0
	maxhh := 0
	
	bounds := Array<Int> new(mesh npolys*4)
	if (!bounds) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMeshDetail: Out of memory 'bounds' (%d).", mesh npolys*4)
		return false
	}
	poly := Array<Float> new(nvp*3)
	if (!poly) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMeshDetail: Out of memory 'poly' (%d).", nvp*3)
		return false
	}
	
	// Find max size for a polygon area.
	for (i: Int in 0..(mesh npolys)) {
		p: UShort* = mesh polys[i*nvp*2]
		// Hmm.. This is where references come in handy (outside of function args)!
		// I really don't wanna do it like this..
		xmin: Int* = bounds[i*4+0]&
		xmax: Int* = bounds[i*4+1]&
		ymin: Int* = bounds[i*4+2]&
		ymax: Int* = bounds[i*4+3]&
		xmin = chf width
		xmax = 0
		ymin = chf height
		ymax = 0
		for (j: Int in 0..nvp) {
			if(p[j] == RC_MESH_NULL_IDX) break
			 v := mesh verts[p[j]*3]&
			// whuuuuuu
			xmin@ = rcMin(xmin, v[0] as Int)
			xmax@ = rcMax(xmax, v[0]as Int)
			ymin@ = rcMin(ymin, v[2] as Int)
			ymax@ = rcMax(ymax, v[2] as Int)
			nPolyVerts+=1
		}
		// *grumble*
		xmin@ = rcMax(0, xmin-1)
		xmax@ = rcMin(chf width, xmax+1)
		ymin@ = rcMax(0, ymin-1)
		ymax@ = rcMin(chf height, ymax+1)
		if (xmin >= xmax || ymin >= ymax) continue
		maxhw = rcMax(maxhw, xmax-xmin)
		maxhh = rcMax(maxhh, ymax-ymin)
	}
	
	hp data = rcAllocArray(UShort, maxhw*maxhh)
	if (!hp data) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMeshDetail: Out of memory 'hp data' (%d).", maxhw*maxhh)
		return false
	}
	
	dmesh nmeshes = mesh npolys
	dmesh nverts = 0
	dmesh ntris = 0
	dmesh meshes = rcAllocArray(UShort, dmesh nmeshes*4)
	if (!dmesh meshes) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMeshDetail: Out of memory 'dmesh.meshes' (%d).", dmesh nmeshes*4)
		return false
	}

	vcap := nPolyVerts+nPolyVerts/2
	tcap := vcap*2

	dmesh nverts = 0
	dmesh verts = rcAllocArray(Float, vcap*3)
	if (!dmesh verts) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMeshDetail: Out of memory 'dmesh.verts' (%d).", vcap*3)
		return false
	}
	dmesh ntris = 0
	dmesh tris = rcAllocArray(UInt8, tcap*4)
	if (!dmesh tris) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMeshDetail: Out of memory 'dmesh.tris' (%d).", tcap*4)
		return false
	}
	
	for (i: Int in 0..(mesh npolys)) {
		p: UShort* = mesh polys[i*nvp*2]&
		
		// Store polygon vertices for processing.
		npoly := 0
		for (j: Int in 0..nvp) {
			if(p[j] == RC_MESH_NULL_IDX) break
			v: UShort* = mesh verts[p[j]*3]&
			poly[j*3+0] = v[0]*cs
			poly[j*3+1] = v[1]*ch
			poly[j*3+2] = v[2]*cs
			npoly+=1
		}
		
		// Get the height data from the area of the polygon.
		hp xmin = bounds[i*4+0]
		hp ymin = bounds[i*4+2]
		hp width = bounds[i*4+1]-bounds[i*4+0]
		hp height = bounds[i*4+3]-bounds[i*4+2]
		getHeightData(chf, p, npoly, mesh verts, hp, stack)
		
		// Build detail mesh.
		nverts := 0
		if (!buildPolyDetail(poly, npoly, sampleDist, sampleMaxError, chf, hp, verts, nverts, tris, edges, samples)) {
			return false
		}
		
		// Move detail verts to world space.
		for (j: Int in 0. nverts) {
			verts[j*3+0] += orig[0]
			verts[j*3+1] += orig[1] + chf ch; // Is this offset necessary?
			verts[j*3+2] += orig[2]
		}
		// Offset poly too, will be used to flag checking.
		for (j: Int in 0..npoly) {
			poly[j*3+0] += orig[0]
			poly[j*3+1] += orig[1]
			poly[j*3+2] += orig[2]
		}
	
		// Store detail submesh.
		ntris := tris size()/4
		
		dmesh meshes[i*4+0] = dmesh nverts
		dmesh meshes[i*4+1] = nverts as UShort
		dmesh meshes[i*4+2] = dmesh ntris
		dmesh meshes[i*4+3] = ntris as UShort
		
		// Store vertices, allocate more memory if necessary.
		if (dmesh nverts+nverts > vcap) {
			while (dmesh nverts+nverts > vcap)
				vcap += 256
				
			newv: Float* = rcAllocArray(Float, vcap*3)
			if (!newv) {
				if (rcGetLog())
					rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMeshDetail: Out of memory 'newv' (%d).", vcap*3)
				return false
			}
			if (dmesh nverts)
				memcpy(newv, dmesh verts, sizeof(Float)*3*dmesh nverts)
			delete [] dmesh verts
			dmesh verts = newv
		}
		for (j: Int in 0. nverts) {
			dmesh verts[dmesh nverts*3+0] = verts[j*3+0]
			dmesh verts[dmesh nverts*3+1] = verts[j*3+1]
			dmesh verts[dmesh nverts*3+2] = verts[j*3+2]
			dmesh nverts+=1
		}
		
		// Store triangles, allocate more memory if necessary.
		if (dmesh ntris+ntris > tcap) {
			while (dmesh ntris+ntris > tcap)
				tcap += 256
			newt := rcAllocArray(UInt8, tcap*4)
			if (!newt) {
				if (rcGetLog())
					rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMeshDetail: Out of memory 'newt' (%d).", tcap*4)
				return false
			}
			if (dmesh ntris)
				memcpy(newt, dmesh tris, sizeof(UInt8)*4*dmesh ntris)
			delete [] dmesh tris
			dmesh tris = newt
		}
		for (j: Int in 0. ntris) {
			t: Int* = tris[j*4]&
			dmesh tris[dmesh ntris*4+0] = t[0] as UInt8
			dmesh tris[dmesh ntris*4+1] = t[1] as UInt8
			dmesh tris[dmesh ntris*4+2] = t[2] as UInt8
			dmesh tris[dmesh ntris*4+3] = getTriFlags(verts[t[0]*3]&, verts[t[1]*3]&, verts[t[2]*3]&, poly, npoly)
			dmesh ntris+=1
		}
	}
	endTime := rcGetPerformanceTimer()
	
	if (rcGetBuildTimes())
		rcGetBuildTimes() buildDetailMesh += rcGetDeltaTimeUsec(startTime, endTime)

	return true
}

rcMergePolyMeshDetails: func(meshes: RCPolyMeshDetail**, nmeshes: Int, mesh: RCPolyMeshDetail@) -> Bool {
	startTime := rcGetPerformanceTimer()
	
	maxVerts := 0
	maxTris := 0
	maxMeshes := 0

	for (i: Int in 0..nmeshes) {
		if (!meshes[i]) continue
		maxVerts += meshes[i] nverts
		maxTris += meshes[i] ntris
		maxMeshes += meshes[i] nmeshes
	}

	mesh nmeshes = 0
	mesh meshes = rcAllocArray(UShort, maxMeshes*4)
	if (!mesh meshes) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMeshDetail: Out of memory 'pmdtl.meshes' (%d).", maxMeshes*4)
		return false
	}

	mesh ntris = 0
	mesh tris = rcAllocArray(UInt8, maxTris*4)
	if (!mesh tris) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMeshDetail: Out of memory 'dmesh.tris' (%d).", maxTris*4)
		return false
	}

	mesh nverts = 0
	mesh verts = rcAllocArray(Float, maxVerts*3)
	if (!mesh verts) {
		if (rcGetLog())
			rcGetLog() log(RCLogCategory RC_LOG_ERROR, "rcBuildPolyMeshDetail: Out of memory 'dmesh.verts' (%d).", maxVerts*3)
		return false
	}
	
	// Merge datas.
	for (i: Int in 0..nmeshes) {
		dm: RCPolyMeshDetail* = meshes[i]
		if (!dm) continue
		for (j: Int in 0..(dm nmeshes)) {
			dst: UShort* = mesh meshes[mesh nmeshes*4]&
			src: UShort* = dm meshes[j*4]&
			dst[0] = mesh nverts+src[0]
			dst[1] = src[1]
			dst[2] = mesh ntris+src[2]
			dst[3] = src[3]
			mesh nmeshes+=1
		}
			
		for (k: Int in 0..(dm nverts)) {
			vcopy(mesh verts[mesh nverts*3]&, dm verts[k*3]&)
			mesh nverts+=1
		}
		for (k: Int in 0..(dm ntris)) {
			mesh tris[mesh ntris*4+0] = dm tris[k*4+0]
			mesh tris[mesh ntris*4+1] = dm tris[k*4+1]
			mesh tris[mesh ntris*4+2] = dm tris[k*4+2]
			mesh tris[mesh ntris*4+3] = dm tris[k*4+3]
			mesh ntris+=1
		}
	}

	endTime := rcGetPerformanceTimer()
	if (rcGetBuildTimes())
		rcGetBuildTimes() mergePolyMeshDetail += rcGetDeltaTimeUsec(startTime, endTime)
	
	return true
}

