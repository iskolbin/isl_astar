[![license](https://img.shields.io/badge/license-public%20domain-blue.svg)]()

isl\_astar
==========
Single header public domain A\* pathfinding algorithm for C. Algorithm implemented
in graph terms, so you need some code to tune it for 2D grids (see example below).
On the other hand, code is straightforward and quite generic, you can use it for
any topology (3D, hexagonal grids, you name it). Implementation uses indirect
binary heap for fast open list operations. 

License
-------
This software is dual-licensed to the public domain and under the following
license: you are granted a perpetual, irrevocable license to copy, modify,
publish, and distribute this file as you see fit.

Author
------
Original A`* is discovered by Peter Hart, Nils Nilsson and Bertram Raphael of
Stanford Research Institute (now SRI International) in 1968. It is an extension
of Edsger Dijkstra's 1959 algorithm.

This library source code is written by Ilya Kolbin, iskolbin@gmail.com in 2016.

original git URL: https://github.com/iskolbin/isl_astar

Limitations
-----------
* Each transition from node to node must yield non-negative cost.
* Number of neighbors is bounded by `ISL_MAX_NEIGHBORS` const (16 by default) 
which can be redefined.

Usage
-----
First include in the one of your source files implemetation macro before
including header file.

```c
#define ISL_ASTAR_IMPLEMENTATION
#include "isl_astar.h"
```

Next define properties struct populated with needed functions. See complete example 
at the end of readme.

```c
isla_properties properties = {your_get_neighbors, your_neighbor_cost, your_estimate_cost};
isla_result result = {0}; // You will need it anyway
```

Neighbor cost function must return cost of transition from current node to selected neighbor.
By default all cost are `double`, which can be redefined by `define ISLA_COST`. Estimate cost
function is the heuristic which returns approximate cost between far nodes. Quality of the 
heuristic directly affects the speed of convergence of the algorithm. For grids it's often a
line distance between cells. 

Also you need to correctly allocate and initialize graph nodes. On the stack it is
as simple as

```c
isla_node nodes[SIZE] = {0};
```

And finally you can actually start using pathfinding.

isla\_find\_path
----------------
Use this function to find path between selected nodes. Prototype is

```c
isla_result isla_find_path( isla_node *start, 
	isla_node *finish, 
	isla_properties *properties, 
	void *userdata )
```

Function returns struct with fields:

`status` - one of the constants:
  * `ISLA_OK` - path was found,
  * `ISLA_BLOCKED` - path is blocked,
  * `ISLA_ERROR_BAD_ALLOC` - error during memory allocation, `malloc` returned `NULL`,
  * `ISLA_ERROR_BAD_REALLOC` - error during memory reallocation, `realloc` returned `NULL`,
  * `ISLA_ERROR_BAD_ARGUMENT` - wrong arguments passed, NULL start or finish or properties;

`path` - if `status == ISLA_OK` then this field will contain vector structure which can be traveresed like:
```c
isla_result result = isla_find_path( start, finish, properties, NULL );
if ( result.status == ISLA_OK ) {
	int i;
	isla_reverse_path( result.path ); // Or you can change for loop
	for ( i = 0; i < result.path->length; i++ ) {
		isla_node *node = result.path->node[i];
		// Do your stuff
		isla_destroy_path( result.path ); // Manually clean path
	}
}
```
Note that resulting path is _reversed_, maybe you should consider using `isla_reverse_path`. Also note
that after successful find you have to release path manually using `isla_destroy_path`. In other cases
when `status != ISLA_OK` all deallocations are made automatically.

isla\_destroy\_path
-------------------

Deallocates path, you need to call it only if correct path was found.


isla\_reverse\_path
-------------------

Because of the implementation resulting path is returned reversed, i.e. `result.path[0] == finish`.
Sometimes it's not desireable, so you can just reverse it. 


isla\_strstatus
---------------

String representation of `isla_status` enum without `ISLA_` prefix, i.e. `isla_status( ISLA_OK )` is `"OK"`.


Example
-------
In this example we implement simple 2D grid, and find path between chosen points.

```c
#define ISL_ASTAR_IMPLEMENTATION 
#include "isl_astar.h"

typedef struct {
	int x;
	int y;
	isla_node *node; // and node->data will cross-reference to containing Cell
} Cell;

typedef struct {
	int width;
	int height;
	Cell *cells;
	const char **level;
} Grid;

Cell *grid_get( Grid *grid, int x, int y ) {
	return (grid->cells + y*grid->width + x);
}

int grid_iswalkable( Grid *grid, int x, int y ) {
	return x >= 0 &&
		x < grid->width &&
		y >= 0 &&
		y < grid->height &&
		grid->level[y][x] != '#'; 
}

void grid_add_ifwalkable( Grid *grid, int x, int y, int *count, isla_node **out ) {
	if ( grid_iswalkable( grid, x, y )) {
		isla_node *node = grid_get( grid, x, y )->node;
		// Note that this check is just small optimization
		if ( node->status != ISLA_NODE_CLOSED ) {
			out[*count] = node;
			*count += 1;
		}
	}
}

// In this example we use rectangular topology with all diagonal neighbors
// Not that it's common when diagonal neighbor is walkable only when at least 1 orthogonal neighbor is
int get_grid_neighbors( isla_node *node, isla_node **out, void *userdata ) {
	Grid *grid = userdata;
	Cell *cell = node->data;	
	int count = 0;
	int x = cell->x;
	int y = cell->y;
	int w = grid->width-1;
	int h = grid->height-1;
	grid_add_ifwalkable( grid, x-1, y, &count, out );
	grid_add_ifwalkable( grid, x+1, y, &count, out );
	grid_add_ifwalkable( grid, x-1, y-1, &count, out );
	grid_add_ifwalkable( grid, x, y-1, &count, out );
	grid_add_ifwalkable( grid, x+1, y-1, &count, out );
	grid_add_ifwalkable( grid, x-1, y+1, &count, out );
	grid_add_ifwalkable( grid, x, y+1, &count, out );
	grid_add_ifwalkable( grid, x+1, y+1, &count, out );
	return count;
}

// For this example we use common euclidean metrics
isla_cost euclidean_cost( isla_node *node1, isla_node *node2, void *data ) {
	Cell *cell1 = node1->data;
	Cell *cell2 = node2->data;
	int dx = cell1->x - cell2->x;
	int dy = cell1->y - cell2->y;
	return sqrt( dx*dx + dy*dy );
}

int main( int argc, char **argv ) {
	const char *level[] = {
		"^.........",
		"......#...",
		"......#...",
		"...####...",
		".........v"
	};
	isla_properties properties = { get_grid_neighbors, euclidean_cost, euclidean_cost };
	const size_t width = 10;
	const size_t height = sizeof( level ) / sizeof( level[0] );
	isla_result result;
	int x0,y0,x1,y1,row,col,i;
	Cell cells[width*height] = {0};
	isla_node nodes[width*height] = {0};
	Grid grid = {width,height,cells,level};

	for ( row = 0; row < height; row++ ) {
		for ( col = 0; col < width; col++ ) {
			int idx = row*width + col;
			cells[idx].x = col;
			cells[idx].y = row;
			cells[idx].node = nodes + idx;
			cells[idx].node->data = cells + idx;
			if (level[row][col] == '^') {
				x0 = col;
				y0 = row;
			} else if ( level[row][col] == 'v' ) {
				x1 = col;
				y1 = row;
			}
		}
	}

	result = isla_find_path( grid_get(&grid,x0,y0)->node, grid_get(&grid,x1,y1)->node, &properties, (void *) &grid); 
	printf( "%s\n", isla_strstatus( result.status ));

	if ( result.status == ISLA_OK ) {
		for ( i = 0; i < result.path->length; i++ ) {
			Cell *cell = result.path->nodes[i]->data;
			printf( "(%d,%d)", cell.x, cell.y );
		}

		// You have to release memory only when anything is OK
		// If path is blocked/memory error accured all cleanups will be run automatically
		isla_destroy_path( result.path );
	}
	return 0;
}
```
