rcAllocArray: func<T>(type: T, count: Int) -> T* {
	return gc_malloc(type size * count)
}

