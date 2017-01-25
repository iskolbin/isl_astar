/* 
 isl_astar - v0.6.0 - public domain library for graph pathfinding
                      no warranty implied; use at your own risk

 author: Ilya Kolbin (iskolbin@gmail.com)
 url: github.com/iskolbin/isl_termbox

 A* algorithm is invented by Peter Hart, Nils Nilsson and Bertram 
 Raphael of Stanford Research Institute (now SRI International) in 1968.
 It is an extension of Edsger Dijkstra's 1959 algorithm.

 LICENSE

 This software is dual-licensed to the public domain and under the following
 license: you are granted a perpetual, irrevocable license to copy, modify,
 publish, and distribute this file as you see fit.
*/

#ifndef ISLA_INCLUDE_ASTAR_H_
#define ISLA_INCLUDE_ASTAR_H_

#include <stddef.h>

#ifdef ISLA_STATIC
	#define ISLA_DEF static
#else
	#define ISLA_DEF extern
#endif

#ifndef ISLA_COST
	typedef double isla_cost;
#else	
	typedef ISLA_COST isla_cost;
#endif
	
#ifndef ISLA_MAX_NEIGHBORS
	#define ISLA_MAX_NEIGHBORS 16
#endif

#if !defined(ISLA_MALLOC)&&!defined(ISLA_REALLOC)&&!defined(ISLA_FREE)
	#include <stdlib.h>
	#define ISLA_MALLOC malloc
	#define ISLA_REALLOC realloc
	#define ISLA_FREE free
#elif !defined(ISLA_MALLOC)||!defined(ISLA_REALLOC)||!defined(ISLA_FREE)
	#error "You must to define ISLA_MALLOC, ISLA_REALLOC, ISLA_FREE to remove stdlib dependency"
#endif

typedef enum {
	ISLA_NODE_DEFAULT,
	ISLA_NODE_OPENED = 1,
	ISLA_NODE_CLOSED = 2,
} isla_node_status;

typedef struct isla_node isla_node;

struct isla_node {
	size_t index;
	isla_cost g;
	isla_cost f;
	isla_node_status status;
	isla_node *parent;
	void *data;
};

typedef isla_node *(*isla_neighbor)( isla_node *, isla_node *, void *userdata );
typedef isla_cost (*isla_cost_fun)( isla_node *, isla_node *, void *userdata );
typedef int (*isla_predicate)( isla_node *, void *userdata );

typedef enum {
	ISLA_OK,
	ISLA_BLOCKED,
	ISLA_ERROR_BAD_ALLOC,
	ISLA_ERROR_BAD_REALLOC,
	ISLA_ERROR_BAD_ARGUMENTS,
} isla_status;

typedef struct {
	isla_node **nodes;
	size_t allocated;
	size_t length;
} isla_path;

typedef struct {
	isla_status status;
	isla_path *path;
} isla_result;

typedef struct {
	isla_neighbor next_neighbor;
	isla_cost_fun eval_cost;
	isla_cost_fun estimate_cost;
	isla_predicate is_finish_node;
	isla_path *cache_used;
	isla_path *cache_open;
} isla_properties;

#ifdef __cplusplus
extern "C" {
#endif

ISLA_DEF isla_result isla_find_path( isla_node *start, isla_node *finish, isla_properties *properties, void *userdata );
ISLA_DEF void isla_destroy_path( isla_path *path );
ISLA_DEF void isla_reverse_path( isla_path *path );
ISLA_DEF const char *isla_strstatus( isla_status status );

#ifdef __cplusplus
}
#endif

#endif // ISLA_INCLUDE_ASTAR_H_

#ifdef ISL_ASTAR_IMPLEMENTATION

// Minimal dynamic vector implementation for path storage
static isla_path *isla__path_create( size_t n ) {
	isla_path *path = ISLA_MALLOC( sizeof *path );
	n = n <= 0 ? 1 : n;
	if ( path != NULL ) {
		path->nodes = ISLA_MALLOC( n * sizeof( *path->nodes ));
		if ( path->nodes != NULL ) {
			path->allocated = n;
			path->length = 0;
		} else {
			free( path );
			path = NULL;
		}
	}
	return path;
}

void isla_destroy_path( isla_path *path ) {
	if ( path != NULL ) {
		ISLA_FREE( path->nodes );
		ISLA_FREE( path );
	}
}

void isla_reverse_path( isla_path *path ) {
	if ( path != NULL ) {
		size_t i;
		size_t half = path->length / 2;
		for ( i = 0; i < half; i++ ) {
			isla_node *tmp = path->nodes[i];
			path->nodes[i] = path->nodes[path->length-i-1];
			path->nodes[path->length-i-1] = tmp;
		}
	}
}

static isla_status isla__path_grow( isla_path *path, size_t newalloc ) {
	isla_node **newnodes = ISLA_REALLOC( path->nodes, sizeof( *path->nodes ) * newalloc );
	if ( newnodes != NULL ) {
		path->allocated = newalloc;
		path->nodes = newnodes;
		return ISLA_OK;
	} else {
		return ISLA_ERROR_BAD_REALLOC;
	}
}

static isla_status isla__path_push( isla_path *path, isla_node *node ) {
	size_t index = path->length;
	isla_status status = ISLA_OK;
	if ( path->allocated <= path->length ) {
		status = isla__path_grow( path, path->allocated * 2 );
		if ( status != ISLA_OK ) {
			return status;
		}
	}
	path->nodes[index] = node;
	path->length++;
	return status;
}
// End of dynamic vector implementation for path storage


// Binary heap implementation for fast open list
static isla_status isla__heap_grow( isla_path *heap, size_t newalloc ) {
	isla_node **newnodes = ISLA_REALLOC( heap->nodes, sizeof( *heap->nodes ) * newalloc );
	if ( newnodes != NULL ) {
		heap->allocated = newalloc;
		heap->nodes = newnodes;
		return ISLA_OK;
	} else {
		return ISLA_ERROR_BAD_ALLOC;
	}
}

static void isla__heap_swap( isla_path *heap, size_t index1, size_t index2 ) {
	isla_node *tmp = heap->nodes[index1];
	heap->nodes[index1] = heap->nodes[index2];
	heap->nodes[index2] = tmp;
	heap->nodes[index1]->index = index1;
	heap->nodes[index2]->index = index2;
}

static size_t isla__heap_siftup( isla_path *heap, size_t index ) {
	size_t parent = (index-1) >> 1;
	while ( index > 0 && heap->nodes[index]->f < heap->nodes[parent]->f ) {
		isla__heap_swap( heap, index, parent );
		index = parent;
		parent = (index-1) >> 1;	
	}
	return index;
}

static void isla__heap_siftdown_floyd( isla_path *heap, size_t index ) {
	size_t left = (index << 1) + 1;
	size_t right = left + 1;
	while ( left < heap->length ) {
		size_t higher = ( right < heap->length && heap->nodes[right]->f < heap->nodes[left]->f ) ? right : left;
		isla__heap_swap( heap, index, higher );
		index = higher;
		left = (index << 1) + 1;
		right = left + 1;
	}
	isla__heap_siftup( heap, index );
}

static void isla__heap_siftdown( isla_path *heap, size_t index ) {
	size_t left = (index << 1) + 1;
	size_t right = left + 1;
	while ( left < heap->length ) {
		size_t higher = ( right < heap->length && heap->nodes[right]->f < heap->nodes[left]->f ) ? right : left;
		if ( heap->nodes[index]->f < heap->nodes[higher]->f ) 
			break;
		isla__heap_swap( heap, index, higher );
		index = higher;
		left = (index << 1) + 1;
		right = left + 1;
	}
}

static isla_status isla__heap_enqueue( isla_path *heap, isla_node *node ) {
	isla_status status = isla__path_push( heap, node );
	if ( status == ISLA_OK ) {
		isla__heap_siftup( heap, heap->length-1 );
	}
	return status;
}

static isla_node *isla__heap_dequeue( isla_path *heap ) {
	if ( heap->length == 0 ) {
		return NULL;
	} else if ( heap->length == 1 ) {
		heap->length = 0;
		return heap->nodes[0];
	} else {
		isla_node *node = heap->nodes[0];
		isla__heap_swap( heap, 0, --heap->length );
		isla__heap_siftdown_floyd( heap, 0 );
		return node;
	}
}

static void isla__heap_update( isla_path *heap, isla_node *node ) {
	isla__heap_siftdown( heap, isla__heap_siftup( heap, node->index ));
}
// End of binary heap implementation


static void isla__cleanup( isla_path *used, isla_path *open, int cached ) {
	size_t i;
	for ( i = 0; i < used->length; i++ ) {
		isla_node *node = used->nodes[i];
		node->g = 0;
		node->f = 0;
		node->status = ISLA_NODE_DEFAULT;
		node->parent = NULL;
		node->index = 0;
	}
	if ( cached ) {
		used->length = 0;
		open->length = 0;
	} else {
		isla_destroy_path( used );
		isla_destroy_path( open );
	}
}

static isla_result isla__build_path( isla_node *node ) {
	isla_result result = {ISLA_OK, isla__path_create(4)};

	if ( result.path == NULL ) {
		result.status = ISLA_ERROR_BAD_ALLOC;
		return result;
	}

	do {
		result.status =  isla__path_push( result.path, node );
		if ( result.status != ISLA_OK ) {
			isla_destroy_path( result.path );
			result.path = NULL;
			return result;
		}
		node = node->parent;
	} while ( node != NULL );

	return result;
}

isla_result isla_find_path( isla_node *start, isla_node *finish, isla_properties *properties, void *userdata ) {
	int cached = 0;
	isla_path *openlist;
	isla_path *usedlist;
	isla_result result = {ISLA_OK,NULL};

	if ( start == NULL || finish == NULL || properties == NULL ) {
		result.status = ISLA_ERROR_BAD_ARGUMENTS;
		return result;
	}

	cached = properties->cache_open != NULL && properties->cache_used != NULL;

	openlist = cached ? properties->cache_open : isla__path_create( ISLA_MAX_NEIGHBORS );
	usedlist = cached ? properties->cache_used : isla__path_create( 4 );

	if ( openlist == NULL || usedlist == NULL ) {
		isla_destroy_path( openlist );
		isla_destroy_path( usedlist );
		result.status = ISLA_ERROR_BAD_ALLOC; 
		return result;
	}

	result.status = isla__heap_enqueue( openlist, start );
	if ( result.status != ISLA_OK ) {
		isla__cleanup( usedlist, openlist, cached );
		return result;
	}

	start->g = 0;
	start->f = properties->estimate_cost( start, finish, userdata );

	while ( openlist->length > 0 ) {
		int i; 
		int count;
	 	isla_node *node = isla__heap_dequeue( openlist );
		isla_node *neighbor = NULL;
		node->status = ISLA_NODE_CLOSED;

		if ( node == finish || (properties->is_finish_node != NULL && properties->is_finish_node( node, userdata ))) {
			result = isla__build_path( finish );
			isla__cleanup( usedlist, openlist, cached );
			return result;
		}

		while ( (neighbor = properties->next_neighbor( node, neighbor, userdata ))) {
			if ( neighbor->status != ISLA_NODE_CLOSED ) {
				isla_cost g = node->g + properties->eval_cost( node, neighbor, userdata );
				if ( neighbor->status == ISLA_NODE_DEFAULT || g < neighbor->g ) {
					neighbor->g = g;
					neighbor->f = g + properties->estimate_cost( node, finish, userdata );
					neighbor->parent = node;
					if ( neighbor->status == ISLA_NODE_OPENED ) {
						isla__heap_update( openlist, neighbor );
					} else {
						neighbor->status = ISLA_NODE_OPENED;
						result.status = isla__heap_enqueue( usedlist, neighbor );
						if ( result.status != ISLA_OK ) {
							isla__cleanup( usedlist, openlist, cached );
							return result;
						}
						result.status = isla__path_push( openlist, neighbor );
						if ( result.status != ISLA_OK ) {
							isla__cleanup( usedlist, openlist, cached );
							return result;
						}
					}
				}
			}
		}
	}

	result.status = ISLA_BLOCKED;
	isla__cleanup( usedlist, openlist, cached );

	return result;
}

const char *isla_strstatus( isla_status status ) {
	const char *statuses[] = {
		"OK",
		"BLOCKED",
		"ERROR_BAD_ALLOC",
		"ERROR_BAD_REALLOC",
		"ERROR_BAD_ARGUMENTS",
	};
	return statuses[status];
}

#endif // ISL_ASTAR_IMPLEMENTATION
