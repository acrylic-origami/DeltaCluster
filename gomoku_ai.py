from functools import reduce
import itertools
from random import random
import time
import math
import sys

colors = {'b':0, 'w':1}
inverse_colors = { 0:'b', 1:'w' }
TIME_LIMIT = 0.25

# the goal: 
    # - efficient expression of board state through disjoint set segments and [hopefully] sparse perimeters that are scored based on density of connections to adjacent segments.
    # - using deltas to represent the game tree rather than copying for deeper searches using less memory and, hopefully, time
        # - the board state is basically global, but I mutate and revert it by travelling up and down the game tree.
    # - the rest is just basic minimax using a BFS on both colors simultaneously. Encountering a trivial trap state will cause backtracking to see if it is non-trivial.

class Delta:
    def __init__(self):
        self.Q = []
        
        self.SET = 0
        self.APPLY = 1
        
        self.inverses = {
            # ("list", "append"): lambda obj: None,
            ("list", "append"): ("pop", ),
            ("list", "pop"): ("append", lambda obj, *args: (obj[-1], )),
            ("dict", "pop"): ("__setitem__", lambda obj, *args: tuple(args + [obj[args[0]]])) #meh, not the most elegant solution, but w/e
        }
    
    def set(self, obj, idx, v):
        gets = [
            lambda x, idx: getattr(x, idx),
            lambda x, idx: x[idx]
        ]
        sets = [
            lambda x, idx, y: setattr(x, idx, y),
            lambda x, idx, y: x.__setitem__(idx, y)
        ] #meh closure
        self.Q.append((self.SET, obj, idx, gets[hasattr(obj, '__getitem__')](obj, idx)))
        sets[hasattr(obj, '__setitem__')](obj, idx, v)
    
    def apply(self, obj, fn, args):
        # grr... we're going to have to build a library of inverse functions >_>
        # e.g. append -> pop
        
        inverse = self.inverses[(type(obj).__name__, fn)]
        if len(inverse)>1:
            #the inverse has arguments: store them
            self.Q.append((self.APPLY, obj, fn, inverse[1](obj, *args)))
        else:
            self.Q.append((self.APPLY, obj, fn))
            
        getattr(obj, fn)(*args) # I guess, instead of using lambdas which would be nice, I have to be explicit with a function name string. MEH.
    
    def merge(self, D):
        self.Q += D.Q
        
    def undo(self):
        for d in self.Q:
            if d[0] == self.SET:
                sets = [
                    lambda x, idx, y: setattr(x, idx, y),
                    lambda x, idx, y: x.__setitem__(idx, y)
                ]
                sets[hasattr(d[1], '__setitem__')](*d[1:])
            # ohhhh boy.
            else:
                fn_inv = getattr(d, self.inverses[(type(d[1]).__name__, d[2])][0])
                if len(d)>4:
                    fn_inv(*d[4])
                else:
                    fn_inv()

def sgn(x):
    return [1, -1][x<=0]
        
def resolve_direction(x, y):
    if y<=0:
        #this applies too to the case (1, -1): we want this to be top-right to bottom-left, not vice versa
        return (-x, -y)
    return (x, y)
        
class DisjointPoint:
    def __init__(self, rep = -1, rank = 0):
        self.rep = rep
        self.rank = rank

class DisjointSet:
    def __init__(self):
        self.d = {}
        
    def find(self, k):
        if self.d[k].rep==k:
            return k
        else:
            self.d[k].rep = self.find(k)
            return self.d[k].rep
    
    def has(self, k):
        return k in self.d
    
    def get(self, idx):
        return self.d[idx]
    
    def union(self, i, j, D = Delta()):
        a = self.d[self.find(i)]
        b = self.d[self.find(j)]
        if a.rank>b.rank:
            D.set(b, 'rep', a)
            try:
                D.merge(a.update_aggregates(self.d[b]))
            except AttributeError:
                pass
                
        else:
            D.set(a, 'rep', b)
            try:
                D.merge(b.update_aggregates(self.d[a]))
            except AttributeError:
                pass
                
            if b.rank == a.rank:
                D.set(b, 'rank', b.rank + 1)
        
        return D
    
    def add(self, location, color, direction, D = Delta()):
        P = Point(location, color, direction, [location, location], (location, direction), 0)
        D.set(self.d, (location, direction), P)
        return D
    # def union(self, i, j):
    #   super(DisjointSet, self).union(i, j)
        #the new parent will have a new shape that we need to use to update the heap for perimeters in those directions
        
        

class Point(DisjointPoint):
    def __init__(self, location, color, direction, bounds = [], *args):
        super().__init__(*args)
        self.location = location
        self.color = color
        self.direction = resolve_direction(*direction)
        self.bounds = bounds
    
    def length(self):
        if self.rep != (self.location + self.direction):
            return self.get(self.find(self.rep)).length()
            
        if not hasattr(self, 'length'):
            self.length = max(abs(self.bounds[0][0]-self.bounds[1][0]), abs(self.bounds[0][1]-self.bounds[1][1]))
        
        return self.length
    
    def update_aggregates(self, x):
        # x is Point
        # unless first element, expect args[0] to contain the most extreme point so far
        # we can look only at the relative position of the leaders relative to each other, instead of looking at individual directions
        if self.location[0]*self.direction[0]<x.location[0]*self.direction[0]:
            #on the "lower" side (left for most, top right on right-to-left diagonal)
            self.bounds[0] = x.bounds[0]
        else:
            self.bounds[1] = x.bounds[1]

class Space:
    def __init__(self, location, adjacencies):
        # (Int)[2], [{Segment}[4]][2]
        self.location = location
        if adjacencies:
            self.adjacencies = ([adjacencies] if type(adjacencies) != list else adjacencies) #black and white adjacencies
        else:
            self.adjacencies = [dict.fromkeys([(1, 0), (0, 1), (1, 1), (-1, 1)])]*2 #initialize the directions, truncate if necessary
        self.score = 0
        
    def set_score(self, score):
        self.score = score
    
    def register_adjacency(self, color, location, obj, D = Delta()):
        # Boolean, (Int)[2], Segment
        direction = resolve_direction([location[i]-self.direction[i] for i in range(2)])
        if direction in self.adjacencies[color]:
            D.apply(self.adjacencies[color][direction], "append", (obj, ))
        return D
        
    def deregister_direction(self, direction, color, D = Delta()):
        # (Int)[2], Boolean
        D.apply(self.adjacencies[color], 'pop', (direction, ))
        return D
        

class Board:
    perimeters = {}
    tiles = [{}, {}]
    flat_tiles = None
    def __init__(self):
        self.perimeters = {} #b, w, dictionary of spaces
        self.tiles = [{}, {}] #b, w, dictionary of tuple coordinates
        self.flat_tiles = None
        self.segments = DisjointSet() #use a disjointset to track the segments that have been recorded already. The segments will have overlap, so the keys need to include a direction components, i.e. (x, y, direction)
        self.promisings = [{},{}] #these are dictionaries of perimeters with scores depending on the length and density of active segments.
        self.sequence_potential = [[], []] # similar to promisings, these are heaps of (potential, dsj_set_key). Potential is defined as the number of moves ot ignorance required if the opponent is to play. Zero and one move of ignorance imply that if you are to play right now, you will win. Having two one-move-ignorance sequences leads to a forced win too, as does three two-move-ignorance, etc. Degeneracies are really just fancy words for only taking into account the combinations of the lowest-value sequences.
        self.parent = self
    
    def initialize_size(self, size):
        #List<Int>[2]
        self.size = size
    
    def resolve(self):
        # used for resolving a board state based on deltas
        return self # if it's on a return basis, use this
        # pass # if it's on a reference basis, use this
            
    def deresolve(self):
        return self
        
    def update_state(self, board):
        #2D List<Int>
        global colors
        # assume that there is only one piece change
        # diff_generator = ((i, j) in self.tiles[pieces[board[i][j]]] for i in range(len(board)) for j in range(len(board)) if board[i][j] != ' ')
        # for l, diff in enumerate(diff_generator):
        #   if diff:
        for i in range(len(board)):
            for j in range(len(board[i])):
                color = colors[board[i][j]]
                if (i, j) not in self.tiles[color]:
                    self.set((i, j), color)
                    self.update_tupletiles()
                    return
    
    def update_tupletiles(self):
        self.tupletiles = tuple(tuple(tile for tile in color) for color in self.tiles)
        # self.tupletiles = tuple(tuple([coord for coord in self.tiles[0][tile]]+([(i, j)] if tile==color else [])) for tile in itertools.chain(*self.tiles)) #yay generator ^_^
    
    def set(self, location, color):
        # Int, Point, Boolean
        
        D = Delta()
        
        self.tiles[color][location] = True
        # it is worth noting that the tiles are added to the disjointset by register_parameters, not here as you might expect. It's more relevant for parameters to be setting direction.
        D.merge(self.register_perimeters(location, color)) #will register regardless of if isolated or not. Returns new perimeter instances
        # IT IS CRUCIAL THIS GOES FIRST: this adds orphaned points to the disjoint set.
        # this also certainly does not add a perimeter where the tile will be placed, so the next loop will no be affected
        if location in self.perimeters:
            # this is one of the perimeters, we can union some shit
            # we also potentially have to update some perimeters, and nullify some opponent sequences
            adjacencies = self.perimeters[location].adjacencies[color]
            for direction in adjacencies:
                #don't assume that for every key set, there is at least one adjacency. If necessary, prune, though to the inconvenience of keys
                if len(adjacency) > 0 and color in [a.color for a in adjacency]:
                    # This is where we used to add sequences; when we were certain of a merge operation. Now, we're adding our piece 
                    # out in the wild open at the top, 
                    # self.segments.add(location, color, direction) #returns Point object; not necessary for Union
                    if len(adjacency) == 2:
                        if adjacency[0].color == adjacency[1].color == color:
                            self.segments.union(adjacencies[direction][0].location, adjacencies[direction][1].location)
                        elif adjacency[1].color == color and adjacency[0].color != color:
                            self.segments.union(adjacencies[direction][1].location, location)
                    
                    # regardless of what happens, merge the 0th element and location if they match; if both sides are the same color, this will get unioned into the larger segment anyways
                    if adjacencies[direction][0] == color:
                        self.segments.union((adjacencies[direction][0].location, direction), (location, direction))
                        root = self.segments.get(self.segments.find(adjacencies[direction]))
                        new_length = root.length()
                        left = (root.bounds[0][0] - root.direction, root.bound[0][1] - root.direction)
                        right = (root.bounds[1][0] + root.direction, root.bound[1][1] + root.direction)
                        
                        heapq.heappush(self.sequence_potential[color], (new_length, left, right))
                        
                        if new_length == 3 and left in self.tiles[color] and right in self.tiles[color] \
                            or newlength == 4 and (left in self.tiles[color]) != (right in self.tiles[color]):
                            #one move of ignorance, or us to play
                            self.almost_traps[color][(location, direction)] = True
                        elif newlength == 4 and left in self.tiles[color] and right in self.tiles[color]:
                            #zero moves of ignorance, irrelevant who plays
                            self.traps[color][(location, direction)] = True
                            
        
        #update externalities
        self.scan(self.nullify_sequences, location, color, 4) #nullify some opponent sequences if necessary
        self.scan(self.update_promisings, location, color, 1) #this updates the promising values of each of the _adjacent_ spaces to this one. This includes blocks of the opposite color, that may now be bounded, as well as blocks of the same color, that are now extended. We can scan the adjacencies of the space before it is deleted in the last step.
        # we may want to update the promisings of both colors, since we're also eventually playing for white
        
        # heapify
        # heapq.heapify(self.promisings[0]) #heapify at the end: it's faster. Promisings is changed in update_promisings, but we don't want to heapify there...
        # heapq.heapify(self.promisings[1])
        
        if location in self.perimeters:
            D.apply(self.perimeters, 'pop', location) #finally let go of the Space
        
        return D
        
    def scan(self, fn, start, source_color, length):
        # Function, (Int)[2], Boolean, Int
        directions = [-1, 0, 1]
        for heading in itertools.permutations(directions, 2):
            #check for at most four tiles in each direction
            fallout = set()
            for i in range(1, length):
                location = (start[0]+i*heading[0], start[1]+i*heading[1])
                if False in [0 <= location[i] < self.size[i] for i in range(2)]:
                    # crossed board boundary
                    break
                
                if location in self.tiles[source_color] or location in self.tiles[not source_color]:
                    # hit a tile: place in fallout
                    fallout.update(self.segments.find((location, ) + (resolve_direction(*heading), ))) # man, tuple concatenation looks so weird...
            
            for f in fallout:
                fn(f, start, source_color) #pass the root node for this function to do its business
    
    def nullify_sequences(self, node, scan_point, source_color, D = Delta()):
        # ((Int[2]), (Int[2])), Point, Boolean
        if node not in self.tiles[not source_color]:
            return # we don't care about our own pieces
        
        bound_dists = ((abs(bound[lr][xy] - scan_point[xy]), lr) for lr, xy in itertools.permutations([0, 1], 2))
        extreme_bound = max(bound_dists)
        direction = resolve_direction(*[sgn(bounds[extreme_bound[1]][i] - scan_point[i]) for i in range(2)])
        minimum_distance = 5 - extreme_bound[0]
        # this is the minimum outer distance on this direction that is needed to still have space to make this sequence
        for i in range(minimum_distance):
            P = (extreme_bound[1][0] + direction[0]*i, extreme_bound[1][1] + direction[1]*i)
            if False in [0 <= P[i] < self.size[i] for i in range(2)] or P in self.tiles[not color]:
                # this space has been bounded
                # by implication, all spaces within this range are null in this direction, but we'll only worry about the ones on the perimeters
                perim_coords = [(P[0]-direction[0], P[1]-direction[1]), (P[0]+direction[0], P[1]+direction[1])]
                for perim_coord in perim_coords:
                    perim = self.perimeters[perim_coord]
                    D.merge(perim.deregister_direction(direction))
                    
                    if len(perim.adjacencies[color][direction]) == 2:
                        lkey = (perim.location[i] - direction[i] for i in range(2)) + direction
                        rkey = (perim.location[i] + direction[i] for i in range(2)) + direction
                        subroot = self.segments.get(self.segments.find(lkey))
                        D.merge(subroot.union(rkey))
                        subroot.length()
                    elif len(self.adjacencies[color][direction]) == 1:
                        perim.adjacencies[color][direction].length()
                    
                    # we can't get zero-length sequences, it's a perimeter
                
                # make the direction for both spaces BARE!
                break
        return D
    
    def update_promisings(self, node, scan_point, source_color):
        # ((Int[2]), (Int[2])), Point, Boolean
        
        raw_direction = [node[0][i] - scan_point[i] for i in range(2)] # since we're scanning with only at most one step in any direction, this will already be normalized
        direction = resolve_direction(*raw_direction)
        root = self.segments.get(self.segments.find(node))
        extermity = [root.bounds[i] + raw_direction[i] for i in range(2)]
        if extermity in perimeters:
            #... #what we do now will depend on how we measure promisings
            if node[0] not in self.tiles[source_color]:
                pass
            else:
                pass
    
    def is_trivial_trap(self, _to_move):
        # if we are to move, then both single zero-ignorance and one-ignorance-moves are traps
        pass
    
    def register_perimeters(self, location, color, D = Delta()):
        directions = [-1, 0, 1]
        for heading in itertools.permutations(range(len(directions)), 2):
            if heading == (0, 0):
                continue
            coord = tuple(location[i]+directions[heading[i]] for i in range(len(heading)))
            direction = tuple(directions[heading[i]] for i in range(len(heading)))
            
            if False in [0 <= coord[i] < self.size[i] for i in range(2)]:
                continue
            
            if not self.segments.has(location + resolve_direction(*direction)):
                D.merge(self.segments.add(location, color, resolve_direction(*direction)))
            
            if coord in self.perimeters:
                # this is a perimeter of something; register the adjacency
                D.merge(self.perimeters[coord].register_adjacency(location, color))
                
                # just to apply the perimeter delta, we have to store the change in the adjacencies of the perimeter as well. Do I have to mock up a whole new slew of Delta classes? Do I have to dirtily use the original classes as delta classes? How can I identify which properties were changed? 
                # I feel like there must be an elegant solution here...
                # okay, well everything is virtually a reference to an object, and most of the time we just need to specify a single object reference (not the reference within the object, but the purest one we can get) and just change that.
                # Trouble is that oftentimes they are stored within objects themselves.
                
                # I could make a generalized class of Deltas, which store a data structure (heap, dict, list, object [reassigning immutable object properties]), an index and a previous value. Any change I wnat to make will be made through this delta object, which records it. We can chain it easily since the changes are simply stored in one big queue. Awesome! Now we've simplified the problem to only four data structures.
                self.segments.get((location, direction)).length()
                
            elif coord not in self.tiles[not color]:
                # this isn't a perimeter of anything; register it and this adjacency
                D.set(self.perimeters, coord, Space(coord, [location]))
        
        return D

class MonotonicBoard(Board):
    def __init__(self, args):
        super().__init__(args)
    
    #methods here will magically make this class only contain the coordinates of a single color
    
class TreeBoard:
    def __init__(self, parent, move, tail_depth = 1, head_depth = 1):
        self.parent = parent #expect Board or TreeBoard object
        self.move = lambda board: board.set(move[0], move[1])
        self.tail_depth = tail_depth
        self.head_depth = head_depth
    
    def resolve_board(self, board):
        self.board = board #exclusively Board object
    
    def set_depths(self, tail, head):
        self.tail_depth = tail
        self.head_depth = head
    
    def apply_deltas(self, board):
        # oh deltas. This is going to be a total pain.
        
        # okay, let's start with a board state delta.
        
        self.D = self.move(board)
        
        #heyyy that was easy!
        return board
    
    def deapply_deltas(self):
        self.D.undo()
        #wow that was also easy... huh.
    
    def resolve(self, d = -1): #automatically resolve from the head down, unless you're part of a recursive chain
        if d != 0:
            self.board = self.apply_deltas(self.parent.resolve(self.head_depth if d == -1 else (d - 1))) #ooooh, the elegance :D
            
            if d == -1: #we might want to do this outside, no?
                self.board.update_tupletiles()
                self.board.update_promisings()
                
            return self.board #this line might not actually be necessary
    
    def deresolve(self, d = -1):
        #Board, Int
        #deresolve to a certain height that is needed to move on to the next board state
        if d != 0:
            self.deapply_deltas(board)
            self.parent.deresolve(board, self.tail_depth if d == -1 else (d - 1))
    
if 'board' not in globals():
    board = Board()

def segmented_BFS(board, color):
    global TIME_LIMIT, board_states, tic
    #return a list of potential moves with a score; we can sort and return
    MAGIC_CONSTANTS = [0.2, 0.05] #we care less about our opponent than we do our moves
    
    # we can cache parts of the tree in a global somehow... that's the utility of running the function even when we're forced to make a move
    
    # it turns out that ancestry isn't needed, since the tree is stored implicitly in TreeBoard anyways.
    
    Q = [board.promisings[color][promising][1] for promising in range(int(len(board.promisings[color])*MAGIC_CONSTANTS[0]))] #filled with Board and [unresolved] TreeBoard states
    Q += [board.promisings[not color][promising][1] for promising in range(int(len(board.promisings[color])*MAGIC_CONSTANTS[1]))]
    level = 1
    while True:
        board_states = {}
        print(board_states)
        tQ = Q[:]
        Q = []
        for tboard in tQ:
            future_board = tboard.resolve() #this also resolves tupletiles [and promisings apparently?]
            if future_board.tupletiles not in board_states:
                future_Q = [promising[1] for promising in future_board.promisings[color][1:math.ceil(len(board)*MAGIC_CONSTANT)]] #this function is currently magic
                future_Q[0].set_depths(tboard.head_depth, 1) #the first element has a long head to apply the deltas, but short tail to undo them
                future_Q[-1].set_depths(1, tboard.tail_depth) #vice versa for the last
                Q += future_Q
                tboard.deresolve() #push to tail
        
        Q[0].set_depths(level, Q[0].tail_depth) #first element must resolve from root
        Q[-1].set_depths(Q[-1].head_depth, level) #last element must undo to root
        
        toc = time.time()
        if (toc-tic) + (toc-tic)/level >= TIME_LIMIT:
            # we aren't going to risk it: break pre-emptively if the average time will make it run over time
            # IF YOU HAVE AN OVERRIDE CONDITION, SPEAK NOW OR FOREVER HOLD YOUR PIECE
            break
        
        level += 1

def get_move(new_state, color):
    global colors, board, tic
    tic = time.time()
    color = colors[color]
    #segmented BFS minimax: we'll start with some fraction of the most promising inner points, then try out some boundaries [based on some other conditions; we can only search to a certain depth before the time is up, so we need to constantly watch the clock]
    board.initialize_size([len(new_state), len(new_state[0])])
    board.update_state(new_state)
    values = segmented_BFS(board, color)
    
    trap_move = board.is_trivial_trap(not color)
    if trap_move:
        return trap_move
        
if __name__ == '__main__':
    state = [
        ['b',' ',' ','b',' ',' ','b',' '],
        [' ','b',' ','b',' ',' ','w',' '],
        [' ','b','b',' ',' ','w','w','w'],
        [' ',' ',' ','w',' ','w','w',' '],
        [' ',' ',' ','w','w','b','b',' '],
        [' ',' ',' ',' ',' ',' ',' ',' '],
        [' ',' ','b','b','w',' ',' ',' '],
        [' ',' ','b','w','w',' ',' ',' ']
    ]
    get_move(state, 'b')