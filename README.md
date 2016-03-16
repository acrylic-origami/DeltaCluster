# &Delta;Cluster: Gomoku AI with clustering heuristic

## Background

[Gomoku](https://en.wikipedia.org/wiki/Gomoku) is most easily described as a generalized version of Tic Tac Toe, whereby both players attempt to create a run of at least _n_ diagonally-, vertically- or horizontally-adjacent pieces on an _m by m_ non-wrapping board. While completely solving Tic Tac Toe is a popular exercise in introductory computational game theory, the naive complexity of a perfect solution increases exponentially with respect to the size of the board, making the development of good intuitive heuristics an appealing strategy for generating a decent AI for general Gomoku with minimal formal study. The advantages of searching a large area within the game tree is, however, not to be understated. With this held constant, the advantage of prioritizing areas of the board to test is visibly vital, so that high potential positions can be better understood by the program, maximizing expected effectiveness. This is in a similar vein to AlphaGo's policy network, which it uses to truncate breadth of its search tree.

## Implementation

This specific AI attempts to fulfill four objectives:

1. Employ a fast clustering heuristic that identifies areas of interest by the quality of the boundaries/liberties (to borrow terminology from Go) of a contiguous block of friendly pieces;
2. Decouple the space complexity and redeuce the copying overhead of tree search from the size of the board by using deltas to represent the game tree traversal rather than the naive solution of passing new board objects to descendants;
3. Specialize in finding complex, multi-move trap configurations;
4. Maximize time consumption within the turn-by-turn time limit.

General AI notes: the AI searches the game tree _without_ making moves for its opponent to find forced wins and trap configurations, then works backwards to determine how likely they are based on chances for the opponent to thwart them by either moving in a space necessary for the trap, or by forcing a move itself. _The value function of the opponent's responses has yet to be implemented._

### Clustering

The speed of the clustering heuristic is owed to the disjoint set structure it uses to test the membership of a new piece to existing cluster[s]. Storing cumulative properties of each cluster becomes fast and easy, making the valuation of clusters and identification of liberties equally so. The latter point is especially important in this version of the algorithm, as the algorithm prunes the sides of game tree by only scanning liberties and the liberties' liberties (for both the AI's color and the opponent's). (currently untested: see `DisjointSet`, `DisjointPoint`, `Space` and `Board` classes)

### Delta

A generalized Delta class is attempted for most common data structure operations in Python, so that board state modifications can be accurately and deterministically undone from the information in Delta objects alone _(currently magic: see `Delta` and `TreeBoard` classes)_

### Trap configurations

Forced wins are propagated up the game tree once they are found, where their "difficulty" (at the current moment, the distance from the last non-forcing move to a forced win) are multiplied by a probability created from value estimations of the player's and opponent's positions to determine the final value of a future position. _This function is currently magic (see `Board::update_promisings`)_

### Time consumption

The program is self-regulating with respect to time. The original competition parameters specified a 250ms time limit per move, and it heeds this by calculating the average time to perform each board analysis, and exiting before the time limit expires assuming this average time cost.

The algorithm is also biased towards evaluating edge-case positions for its own pieces, searching 20% of any node's children, in contrast to the opponent position to which it only searched 5%. The algorithm was designed primarily to best naive single-move-lookahead AIs that played at worst one-move-to-win trap configurations.