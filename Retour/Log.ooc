RCLogCategory: class {
	RC_LOG_PROGRESS: static const Int = 1
	RC_LOG_WARNING: static const Int = 2
	RC_LOG_ERROR: static const Int = 3
}

RCLog: class {
	init: func() {}
	
	log: func(category: Int, format: String, ...) {
		if (m_messageCount >= MAX_MESSAGES)
			return
		
		list: VaList
		va_start(list, format)
        length := vsnprintf(null, 0, format, list) + 1
        output: String = gc_malloc(length)
        va_end(list)
		
        va_start(list, this)
        vsnprintf(output, length, format, list)
        va_end(list)
		m_messages[m_messageCount] = RCLogMessage new(category, output)
		m_messageCount += 1
	}
	
	clear: inline func() {
		m_messageCount = 0
	}
	
	getMessageCount: inline func() -> Int {
		return m_messageCount
	}
	getMessageType: inline func(i: Int) -> Int {
		return m_messages[i] type
	}
	getMessageText: inline func(i: Int) -> String {
		return m_messages[i] message
	}
	
	MAX_MESSAGES: static const Int = 1000
	m_messages: RCLogMessage[RcLog MAX_MESSAGES]
	m_messageCount: Int
}

RCLogMessage: class {
	type: Int
	message: String
	init: func(=type, =message) {}
}

RCBuildTimes: class {
	rasterizeTriangles: Int
	buildCompact: Int
	buildContours: Int
	buildContoursTrace: Int
	buildContoursSimplify: Int
	filterBorder: Int
	filterWalkable: Int
	filterMarkReachable: Int
	buildPolymesh: Int
	erodeArea: Int
	buildDistanceField: Int
	buildDistanceFieldDist: Int
	buildDistanceFieldBlur: Int
	buildRegions: Int
	buildRegionsReg: Int
	buildRegionsExp: Int
	buildRegionsFlood: Int
	buildRegionsFilter: Int
	buildDetailMesh: Int
	mergePolyMesh: Int
	mergePolyMeshDetail: Int
}

g_log: static RCLog = null;
g_btimes: static RCBuildTimes =  null;

rcSetLog: func(log: RCLog) {
	g_log = log
}

rcGetLog: func() -> RCLog {
	return g_log
}

rcSetBuildTimes: func(btimes: RCBuildTimes) {
	g_btimes = btimes
}

rcGetBuildTimes: func() -> RCBuildTimes {
	return g_btimes
}

