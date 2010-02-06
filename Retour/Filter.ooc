use retour
import Retour/[Retour, Timer]

// TODO: Missuses ledge flag, must be called before rcFilterLedgeSpans!
rcFilterLowHangingWalkableObstacles: func(walkableClimb: Int, solid: RCHeightfield@) {
	w := solid width
	h := solid height
	
	for (y: Int in 0..h) {
		for (x: Int in 0..w) {
			ps: RCSpan* = 0
			for (s: RCSpan* = solid spans[x + y*w]; s; s = (s next)) {
				walkable: Bool = (s flags & RC_WALKABLE) != 0
				previousWalkable: Bool = ps && (ps flags & RC_WALKABLE) != 0
				// If current span is not walkable, but there is walkable
				// span just below it, mark the span above it walkable too.
				// Missuse the edge flag so that walkable flag cannot propagate
				// past multiple non-walkable objects.
				if (!walkable && previousWalkable) {
					if (rcAbs((s smax as Int) - (ps smax as Int)) <= walkableClimb)
						s flags |= RC_LEDGE
				}
				ps = s
			}
			// Transfer "fake ledges" to walkables.
			for (s: RCSpan* = solid spans[x + y*w]; s; s = s next) {
				if (s flags & RC_LEDGE)
					s flags |= RC_WALKABLE
				s flags &= ~RC_LEDGE
				ps = s
			}
		}
	}
}
	
rcFilterLedgeSpans: func(walkableHeight: Int, alkableClimb: Int, solid: RCHeightfield@) {
	startTime := rcGetPerformanceTimer()
	
	w := solid width
	h := solid height
	MAX_HEIGHT: static const Int = 0xffff
	
	// Mark border spans.
	for (y: Int in 0..h) {
		for (x: Int in 0..w) {
			for (s: RCSpan* = solid spans[x + y*w]; s; s = s next) {
				// Skip non walkable spans.
				if ((s flags & RC_WALKABLE) == 0)
					continue
				
				bot := (s smax as Int)
				top := s next ? (s next smin as Int) : MAX_HEIGHT
				
				// Find neighbours minimum height.
				minh := MAX_HEIGHT

				// Min and max height of accessible neighbours.
				asmin := s smax
				asmax := s smax
				
				for (dir: Int in 0..4) {
					dx := x + rcGetDirOffsetX(dir)
					dy := y + rcGetDirOffsetY(dir)
					// Skip neighbours which are out of bounds.
					if (dx < 0 || dy < 0 || dx >= w || dy >= h) {
						minh = rcMin(minh, (-walkableClimb) - bot)
						continue
					}
					
					// From minus infinity to the first span.
					ns: RCSpan* = solid spans[dx + dy*w]
					nbot := -walkableClimb
					ntop := ns ? (ns smin as Int) : MAX_HEIGHT
					// Skip neightbour if the gap between the spans is too small.
					if (rcMin(top,ntop) - rcMax(bot,nbot) > walkableHeight)
						minh = rcMin(minh, nbot - bot)
					
					// Rest of the spans.
					for (ns = solid spans[dx + dy*w]; ns; ns = ns next) {
						nbot = ns smax as Int
						ntop = ns next ? (ns next smin as Int) : MAX_HEIGHT
						// Skip neightbour if the gap between the spans is too small.
						if (rcMin(top, ntop) - rcMax(bot, nbot) > walkableHeight) {
							minh = rcMin(minh, nbot - bot)
							// Find min/max accessible neighbour height. 
							if (rcAbs(nbot - bot) <= walkableClimb) {
								if (nbot < asmin) asmin = nbot
								if (nbot > asmax) asmax = nbot
							}
							
						}
					}
				}
				
				// The current span is close to a ledge if the drop to any
				// neighbour span is less than the walkableClimb.
				if (minh < -walkableClimb)
					s flags |= RC_LEDGE
					
				// If the difference between all neighbours is too large,
				// we are at steep slope, mark the span as ledge.
				if ((asmax - asmin) > walkableClimb)
					s flags |= RC_LEDGE
			}
		}
	}
	
	endTime := rcGetPerformanceTimer()
//	if (rcGetLog())
//		rcGetLog()->log(RC_LOG_PROGRESS, "Filter border: %.3f ms", rcGetDeltaTimeUsec(startTime, endTime)/1000.0f)
	if (rcGetBuildTimes())
		rcGetBuildTimes() filterBorder += rcGetDeltaTimeUsec(startTime, endTime)
}	

rcFilterWalkableLowHeightSpans: func(walkableHeight: Int, solid: RCHeightfield@) {
	startTime := rcGetPerformanceTimer()
	
	w := solid width
	h := solid height
	MAX_HEIGHT: static const Int = 0xffff
	
	// Remove walkable flag from spans which do not have enough
	// space above them for the agent to stand there.
	for (y: Int in 0..h) {
		for (x: Int in 0..w) {
			for (s: RCSpan* = solid spans[x + y*w]; s; s = s next) {
				bot = (s smax as Int)
				top = s next ? (s next smin as Int) : MAX_HEIGHT
				if ((top - bot) <= walkableHeight)
					s flags &= ~RC_WALKABLE
			}
		}
	}
	
	endTime := rcGetPerformanceTimer()

//	if (rcGetLog())
//		rcGetLog()->log(RC_LOG_PROGRESS, "Filter walkable: %.3f ms", rcGetDeltaTimeUsec(startTime, endTime)/1000.0f)
	if (rcGetBuildTimes())
		rcGetBuildTimes() filterWalkable += rcGetDeltaTimeUsec(startTime, endTime)
}

