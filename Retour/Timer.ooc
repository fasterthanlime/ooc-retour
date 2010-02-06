rcGetPerformanceTimer: func() -> LLong {
	return Time microtime()
}

rcGetDeltaTimeUsec: func(start, end: LLong) -> Int {
	return (end - start) as Int
}

