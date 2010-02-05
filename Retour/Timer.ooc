rcGetPerformanceTimer: func() -> LLong {
	return Time microtime()
}

rcGetDeltaTimeUsec(start, end: LLong) -> Int {
	return (end - start) as Int
}

