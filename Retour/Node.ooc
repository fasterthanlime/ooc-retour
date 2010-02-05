
static const UShort DT_NULL_IDX = 0xffff

//enum
DTNodeFlags: class {
	DT_NODE_OPEN: static const Int = 0x01
	DT_NODE_CLOSED: static const Int = 0x02
}

DTNode: class {
	cost: Float
	total: Float
	id: UInt
	pidx: UInt //: 30
	flags: UInt //: 2
}

DTNodePool: class {
	DTNodePool(=m_maxNodes, =m_hashSize) {
		m_nodes = rcAllocArray(DTNode, m_maxNodes)
		m_next = rcAllocArray(UShort, m_maxNodes)
		m_first = rcAllocArray(UShort, m_hashSize)
		memset(m_first, 0xff, sizeof(UShort)*m_hashSize)
		memset(m_next, 0xff, sizeof(UShort)*m_maxNodes)
	}

	destroy: func() {
		delete [] m_nodes
		delete [] m_next
		delete [] m_first
	}

	clear: func() {
		memset(m_first, 0xff, sizeof(UShort)*m_hashSize)
		m_nodeCount = 0
	}

	findNode: func(id: UInt) -> DTNode* {
		bucket: UInt = hashint(id) & (m_hashSize-1)
		i := m_first[bucket]
		while (i != DT_NULL_IDX) {
			if (m_nodes[i] id == id)
				return &m_nodes[i]
			i = m_next[i]
		}
		return 0
	}

	getNode: func(id: UInt) -> DTNode* {
		bucket: UInt = hashint(id) & (m_hashSize-1)
		i := m_first[bucket]
		node: DTNode* = 0
		while (i != DT_NULL_IDX) {
			if (m_nodes[i] id == id)
				return &m_nodes[i]
			i = m_next[i]
		}
		
		if (m_nodeCount >= m_maxNodes)
			return 0
		
		i = m_nodeCount as UShort
		m_nodeCount++
		
		// Init node
		node = &m_nodes[i]
		node pidx = 0
		node cost = 0
		node total = 0
		node id = id
		node flags = 0
		
		m_next[i] = m_first[bucket]
		m_first[bucket] = i
		return node
	}
	
	getNodeIdx: inline func(node: DTNode*) -> UInt {
		if (!node) return 0
		return ((node - m_nodes)+1 as UInt)
	}
	
	getNodeAtIdx: inline func(idx: UInt) -> DTNode* {
		if (!idx) return 0
		return &m_nodes[idx-1]
	}
	
	getMemUsed: inline func() -> Int {
		return sizeof(This) +
		sizeof(DTNode)*m_maxNodes +
		sizeof(UShort)*m_maxNodes +
		sizeof(UShort)*m_hashSize
	}
	
//private
	hashint: inline func(UInt a) -> UInt {
		a += ~(a<<15)
		a ^=  (a>>10)
		a +=  (a<<3)
		a ^=  (a>>6)
		a += ~(a<<11)
		a ^=  (a>>16)
		return a
	}
	
	m_nodes: DTNode*
	m_first: UShort*
	m_next: UShort*
	m_maxNodes: Int
	m_hashSize: Int
	m_nodeCount: Int
}

DTNodeQueue: class {
	init: func(=m_capacity) {
		m_heap = rcAllocArray(DTNode*, m_capacity+1)
	}
	
	destroy: func() {
		delete [] m_heap
	}
	
	clear: inline func() {
		m_size = 0
	}
	
	top: inline func() -> DTNode* {
		return m_heap[0]
	}
	
	pop: inline func() -> DTNode* {
		result: DTNode* = m_heap[0]
		m_size--
		trickleDown(0, m_heap[m_size])
		return result
	}
	
	push: inline func(node: DTNode*) {
		m_size++
		bubbleUp(m_size-1, node)
	}
	
	modify: inline func(node: DTNode*) {
		for (i: Int in 0..m_size) {
			if (m_heap[i] == node) {
				bubbleUp(i, node)
				return
			}
		}
	}
	
	empty: inline func() -> Bool {
		return m_size == 0
	}
	
	getMemUsed: inline func() -> Int {
		return sizeof(This) +
		sizeof(DTNode*)*(m_capacity+1)
	}
	
//private
	bubbleUp: func(Int i, node: DTNode*) {
		parent: Int = (i-1)/2
		// note: (index > 0) means there is a parent
		while ((i > 0) && (m_heap[parent] total > node total)) {
			m_heap[i] = m_heap[parent]
			i = parent
			parent = (i-1)/2
		}
		m_heap[i] = node
	}
	
	trickleDown: func(Int i, node: DTNode*) {
		child: Int = (i*2)+1
		while (child < m_size) {
			if (((child+1) < m_size) && 
				(m_heap[child] total > m_heap[child+1] total)) {
				child++
			}
			m_heap[i] = m_heap[child]
			i = child
			child = (i*2)+1
		}
		bubbleUp(i, node)
	}
	
	m_heap: DTNode**
	m_capacity, m_size: Int
}

