
RCConstants: class {
	PI: static const Double = 3.1415926535897932384626
}

rcAllocArray: func<T>(type: T, count: Int) -> T* {
	return gc_malloc(T size * count)
}

// Common helper functions
rcSwap: inline func<T>(a, b: T@) { t: T = a; a = b; b = t }
rcMin: inline func<T>(a, b: T) -> T { return a < b ? a : b }
rcMax: inline func<T>(a, b: T) -> T { return a > b ? a : b }
rcAbs: inline func<T>(a: T) -> T { return a < 0 ? -a : a }
rcSqr: inline func<T>(a: T) -> T { return a * a }
rcClamp: inline func<T>(v, mn, mx: T) -> T { return v < mn ? mn : (v > mx ? mx : v) }

// Common vector helper functions.
vcross: inline func(dest, v1, v2: Float*) {
	dest[0] = v1[1]*v2[2] - v1[2]*v2[1]
	dest[1] = v1[2]*v2[0] - v1[0]*v2[2]
	dest[2] = v1[0]*v2[1] - v1[1]*v2[0]
}

vdot: inline func(v1, v2: Float*) -> Float {
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
}

vmad: inline func(dest, v1, v2: Float*, s: Float) {
	dest[0] = v1[0]+v2[0]*s
	dest[1] = v1[1]+v2[1]*s
	dest[2] = v1[2]+v2[2]*s
}

vadd: inline func(dest, v1, v2: Float*) {
	dest[0] = v1[0]+v2[0]
	dest[1] = v1[1]+v2[1]
	dest[2] = v1[2]+v2[2]
}

vsub: inline func(dest, v1, v2: Float*) {
	dest[0] = v1[0]-v2[0]
	dest[1] = v1[1]-v2[1]
	dest[2] = v1[2]-v2[2]
}

vmin: inline func(mn, v: Float*) {
	mn[0] = rcMin(mn[0], v[0])
	mn[1] = rcMin(mn[1], v[1])
	mn[2] = rcMin(mn[2], v[2])
}

vmax: inline func(mx, v: Float*) {
	mx[0] = rcMax(mx[0], v[0])
	mx[1] = rcMax(mx[1], v[1])
	mx[2] = rcMax(mx[2], v[2])
}

vcopy: inline func(dest, v: Float*) {
	dest[0] = v[0]
	dest[1] = v[1]
	dest[2] = v[2]
}

vdist: inline func(v1, v2: Float*) -> Float {
	dx := v2[0] - v1[0]
	dy := v2[1] - v1[1]
	dz := v2[2] - v1[2]
	return sqrt(dx*dx + dy*dy + dz*dz)
}

vdistSqr: inline func(v1, v2: Float*) -> Float {
	dx := v2[0] - v1[0]
	dy := v2[1] - v1[1]
	dz := v2[2] - v1[2]
	return dx*dx + dy*dy + dz*dz
}

vnormalize: inline func(v: Float*) {
	d: Float = 1.0 / sqrt(rcSqr(v[0]) + rcSqr(v[1]) + rcSqr(v[2]))
	v[0] *= d
	v[1] *= d
	v[2] *= d
}

vequal: inline func(p0, p1: Float*) -> Bool {
	thr: static Float = rcSqr(1.0/16384.0)
	d := vdistSqr(p0, p1)
	return d < thr
}

nextPow2: inline func(v: UInt) -> UInt {
	v -= 1
	v |= v >> 1
	v |= v >> 2
	v |= v >> 4
	v |= v >> 8
	v |= v >> 16
	v += 1
	return v
}

ilog2: inline func(v: UInt) -> UInt {
	r, shift: UInt
	r = (v > 0xffff) << 4; v = v >> r	// >>= isn't working :'(
	shift = (v > 0xff) << 3; v = v >> shift; r |= shift
	shift = (v > 0xf) << 2; v = v >> shift; r |= shift
	shift = (v > 0x3) << 1; v = v >> shift; r |= shift
	r |= (v >> 1)
	return r
}

align4: inline func(x: Int) -> Int {
	return (x+3) & ~3
}

vdot2D: inline func(u, v: Float*)  -> Float {
	return u[0]*v[0] + u[2]*v[2]
}

vperp2D: inline func(u, v: Float*) -> Float {
	return u[2]*v[0] - u[0]*v[2]
}

triArea2D: inline func(a, b, c: Float*) -> Float {
	return ((b[0]*a[2] - a[0]*b[2]) + (c[0]*b[2] - b[0]*c[2]) + (a[0]*c[2] - c[0]*a[2])) * 0.5
}

checkOverlapBox: inline func(amin, amax, bmin, bmax: UShort[3]) -> Bool {
	overlap := true
	overlap = (amin[0] > bmax[0] || amax[0] < bmin[0]) ? false : overlap
	overlap = (amin[1] > bmax[1] || amax[1] < bmin[1]) ? false : overlap
	overlap = (amin[2] > bmax[2] || amax[2] < bmin[2]) ? false : overlap
	return overlap
}

overlapBounds: inline func(amin, amax,  bmin, bmax: Float*) -> Bool {
	overlap := true
	overlap = (amin[0] > bmax[0] || amax[0] < bmin[0]) ? false : overlap
	overlap = (amin[1] > bmax[1] || amax[1] < bmin[1]) ? false : overlap
	overlap = (amin[2] > bmax[2] || amax[2] < bmin[2]) ? false : overlap
	return overlap
}

closestPtPointTriangle: func(closest, p, a, b, c: Float*) {
	// Check if P in vertex region outside A
	ab, ac, ap: Float[3]
	vsub(ab, b, a)
	vsub(ac, c, a)
	vsub(ap, p, a)
	d1: Float = vdot(ab, ap)
	d2: Float = vdot(ac, ap)
	if (d1 <= 0.0 && d2 <= 0.0) {
		// barycentric coordinates (1,0,0)
		vcopy(closest, a)
		return
	}
	
	// Check if P in vertex region outside B
	bp: Float[3]
	vsub(bp, p, b)
	d3: Float = vdot(ab, bp)
	d4: Float = vdot(ac, bp)
	if (d3 >= 0.0 && d4 <= d3) {
		// barycentric coordinates (0,1,0)
		vcopy(closest, b)
		return
	}
	
	// Check if P in edge region of AB, if so return projection of P onto AB
	vc: Float = d1*d4 - d3*d2
	if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
		// barycentric coordinates (1-v,v,0)
		v: Float = d1 / (d1 - d3)
		closest[0] = a[0] + v * ab[0]
		closest[1] = a[1] + v * ab[1]
		closest[2] = a[2] + v * ab[2]
		return
	}
	
	// Check if P in vertex region outside C
	cp: Float[3]
	vsub(cp, p, c)
	d5: Float = vdot(ab, cp)
	d6: Float = vdot(ac, cp)
	if (d6 >= 0.0 && d5 <= d6) {
		// barycentric coordinates (0,0,1)
		vcopy(closest, c)
		return
	}
	
	// Check if P in edge region of AC, if so return projection of P onto AC
	vb: Float = d5*d2 - d1*d6
	if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
		// barycentric coordinates (1-w,0,w)
		w: Float = d2 / (d2 - d6)
		closest[0] = a[0] + w * ac[0]
		closest[1] = a[1] + w * ac[1]
		closest[2] = a[2] + w * ac[2]
		return
	}
	
	// Check if P in edge region of BC, if so return projection of P onto BC
	va: Float = d3*d6 - d5*d4
	if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
		// barycentric coordinates (0,1-w,w)
		w: Float = (d4 - d3) / ((d4 - d3) + (d5 - d6))
		closest[0] = b[0] + w * (c[0] - b[0])
		closest[1] = b[1] + w * (c[1] - b[1])
		closest[2] = b[2] + w * (c[2] - b[2])
		return
	}
	
	// P inside face region. Compute Q through its barycentric coordinates (u,v,w)
	denom: Float = 1.0 / (va + vb + vc)
	v: Float = vb * denom
	w: Float = vc * denom
	closest[0] = a[0] + ab[0] * v + ac[0] * w
	closest[1] = a[1] + ab[1] * v + ac[1] * w
	closest[2] = a[2] + ab[2] * v + ac[2] * w
}

intersectSegmentPoly2D: func(p0, p1, verts: Float*, nverts: Int, tmin, tmax: Float@, segMin, segMax: Int@) -> Bool {
	EPS: static const Float = 0.00000001
	
	tmin = 0
	tmax = 1
	segMin = -1
	segMax = -1
	
	dir: Float[3]
	vsub(dir, p1, p0)
	
	j := nverts-1
	for (i: Int in 0..nverts) {
		edge, diff: Float[3]
		vsub(edge, verts[i*3]&, verts[j*3]&)
		vsub(diff, p0, verts[j*3]&)
		n := vperp2D(edge, diff)
		d := -vperp2D(edge, dir)
		if (rcAbs(d) < EPS) {
			// S is nearly parallel to this edge
			if (n < 0)
				return false
			else
				continue
		}
		t: Float = n / d
		if (d < 0) {
			// segment S is entering across this edge
			if (t > tmin) {
				tmin = t
				segMin = j
				// S enters after leaving polygon
				if (tmin > tmax)
					return false
			}
		} else {
			// segment S is leaving across this edge
			if (t < tmax) {
				tmax = t
				segMax = j
				// S leaves before entering polygon
				if (tmax < tmin)
					return false
			}
		}
		j = i
	}
	return true
}

distancePtSegSqr2D: func(pt, p, q: Float*, t: Float@) -> Float {
	pqx: Float = q[0] - p[0]
	pqz: Float = q[2] - p[2]
	dx: Float = pt[0] - p[0]
	dz: Float = pt[2] - p[2]
	d: Float = pqx*pqx + pqz*pqz
	t = pqx*dx + pqz*dz
	if (d > 0) t /= d
	if (t < 0) t = 0
	else if (t > 1) t = 1
	dx = p[0] + t*pqx - pt[0]
	dz = p[2] + t*pqz - pt[2]
	return dx*dx + dz*dz
}

calcPolyCenter: func(tc: Float*, idx: UShort*, nidx: Int, verts: Float*) {
	tc[0] = 0.0
	tc[1] = 0.0
	tc[2] = 0.0
	for (j: Int in 0..nidx) {
		v: Float* = verts[idx[j]*3]&
		tc[0] += v[0]
		tc[1] += v[1]
		tc[2] += v[2]
	}
	s: Float = 1.0 / nidx
	tc[0] *= s
	tc[1] *= s
	tc[2] *= s
}

closestHeightPointTriangle: func(p, a, b, c: Float*, h: Float@) -> Bool {
	v0, v1, v2: Float[3]
	vsub(v0, c, a)
	vsub(v1, b, a)
	vsub(v2, p, a)
	
	dot00 := vdot2D(v0, v0)
	dot01 := vdot2D(v0, v1)
	dot02 := vdot2D(v0, v2)
	dot11 := vdot2D(v1, v1)
	dot12 := vdot2D(v1, v2)
	
	// Compute barycentric coordinates
	invDenom: Float = 1.0 / (dot00 * dot11 - dot01 * dot01)
	u: Float = (dot11 * dot02 - dot01 * dot12) * invDenom
	v: Float = (dot00 * dot12 - dot01 * dot02) * invDenom

	// The (sloppy) epsilon is needed to allow to get height of points which
	// are interpolated along the edges of the triangles.
	EPS: static const Float = 1e-4
	
	// If point lies inside the triangle, return interpolated ycoord.
	if (u >= -EPS && v >= -EPS && (u+v) <= 1+EPS) {
		h = a[1] + v0[1]*u + v1[1]*v
		return true
	}
	return false
}

distancePtPolyEdgesSqr: func(pt, verts: Float*, nverts: Int, ed, et: Float*) -> Bool {
	// TODO: Replace pnpoly with triArea2D tests?
	c := false
	j := nverts-1
	for (i: Int in 0..nverts) {
		vi: Float* = verts[i*3]&
		vj: Float* = verts[j*3]&
		if (((vi[2] > pt[2]) != (vj[2] > pt[2])) &&
			(pt[0] < (vj[0]-vi[0]) * (pt[2]-vi[2]) / (vj[2]-vi[2]) + vi[0]))
			c = !c
		ed[j] = distancePtSegSqr2D(pt, vj, vi, et[j])
		j = i
	}
	return c
}

