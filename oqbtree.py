import numpy as np
"""A Hierarchical Adaptive Meshmesh indexing tranlation and rotation space.

Combining octree (for indexing tranlation space) and quad-tree (for indexing sphere) to indexing all configurations of two rigid body.
        ^
       Z|1--+
        |_______
       /|      /|3-++
      / | +++ / |
     /__|___7/  |
 +-+5|  |____|__|__Y__>
     | /0--- | /2-+-
     |/      |/
     |_______|
   X/4+--    6++-
   v       
   

"""


# This dictionary is used by the findBranch function, to return the correct branch index
DIRLOOKUP = {'+++':7, '-++':3, '--+':1, '+-+':5, '++-':6, '-+-':2, '---':0, '+--':4}
#### End Globals ####

class conf:
    def __init__(self, node_idx, position=None, vector=None, angle = None,value=0.0):
        self.idx = node_idx
        self.position = position
        self.vector = vector
        self.angle = angle
        self.value = value

class Node:
    """Node.

    """
    def __init__(self,node_idx, centre, size, leaf_num):
        self.error = 100.0
        self.centre = centre
        self.size = size
        self.isLeafNode = True
        self.idx = node_idx
        self.data = data
        self.leaf_num = leaf_num
        # children store region
        self.children = [None for i in range(self.leaf_num)] # the branches should have order
        # grids store mesh grids
        self.grids = []
        # neighbor store neighbor for iterpolation  in this region.
        self.dict = {}

class Bitree:
    """Bittree

    """
    def __init__(self,node_idx, centre, size=np.pi * 2, position, direction):
        self.position = position
        self.direct = direction
        self.root = Node(node_idx, centre, size, leaf_num=2)
        self.root.grids.append(conf(node_idx+'C1',self.position, self.direct, self.centre - size/2.0, 0.0))
        self.root.grids.append(conf(node_idx+'C0',self.position, self.direct, self.centre, 0.0))
        self.root.grids.append(conf(node_idx+'C1',self.position, self.direct, self.centre + size/2, 0.0))
        self.nodes = {}
        self.grids = {}
        self.iterateGrid()
        self.iterateNode()
    def subdivideNode(self, parent):
        parent.isLeafNode = False
        _offset = parent.size/4.0
        _centre = (parent.centre - _offset, parent.centre + _offset )
        for i in range(2):
            parent.children.append(Node(parent.idx+str(i),_centre[i], self.size/2.0, 2))
        left = parent.children[0]
        left.grids.append(parent.grids[0])
        left.grids.append(conf(left.idx + 'C0', self.position, self.direct, left.centre, 0.0))
        left.grids.append(parent.grids[1])
        right = parent.children[1]
        right.grids.append(parent.grids[1])
        right.grids.append(conf(left.idx + 'C0', self.position, self.direct, right.centre, 0.0))
        right.grids.append(parent.grids[2])
	self.iterateGrid()
	self.iterateNode()

    def fill(self, conf_idx, value):
        """fill conf after generation 
        """
        if conf_idx in self.grids:
            self.grids[conf_idx].value = value
        else:
            self.addNode(conf_idx)
            if conf_idx in self.grids:
                self.grids[conf_idx].value = value
            else:
                raise Exception("Con't fill conf %s\n"%(conf_idx))

    def addNode(self,node_idx):
        node_idx = node_idx.split('C')[0]
        pre_idx, idxs = node_idx.split('N')
        pre_idx += 'N' 
        for i in idxs:
            pre_idx +=  str(i)
            if pre_idx not in self.nodes:
                self.subdivideNode(self.nodes[pre_idx[:-1]])
                self.iterateNode()

    def interpolation(self, angle):
        neighbors = self.findNeighbors(self.root, angle)
        v1 = neighbors[0].value
        v2 = neighbors[1].value
        w1 = angle - neighbors[0].angle
        w2 = neighbors[1].angle -angle
        value = (v1 * w1 + v2 * w2)/(w1 + w2)
        return value

    def findNeighbors(self, node, angle):
        _neighbor = None
        if node == None:
            return None  
        elif node.isLeafNode:
            _neighbor = (node.grids[0], node.grids[2])
            return _neighbor
        else:
            child = self.findChild(node, angle)
            return findNeighbors(node.children[child],angle)
            
    def findChild(self, node, angle):
        child_idx = None
        if angle < node.centre: 
            child_idx = 0
        else:
            child_idx = 1
        return  child_idx

    def iterateGrid(self):
        for conf in self._iterateGrid_help(self.root):
            if conf.idx not in self.grids:
                self.grids[conf.idx] = conf
        
    def _iterateGrid_help(self, node):
        """iterate all conf, not unique

        """
        for conf in node.grids:
            if conf != None:yield conf
        for child in node.children:
            for c in self._iterateGrid_help(child):
                yield c 

    def iterateNode(self):
        for n in self._iterateNode_hel p(self.root):
            if n.idx not in self.nodes:
                self.nodes[n.idx] = n

    def _iterateNode_help(self, node):
        if node != None:
            yield node 
        for child in node.children:
            for n in self._iterateNode_help(child):
                yield  n
            
class Quadtree:
    """QuadTree to index sphere

    """
    def __init__(self, com, node_idx):
        """init 6 grids and 8 triangles in root node.
        
        com is centre of mass.
        grid order is like '4', triangles child is anti o' clock from (1,1,1)
        """
        self.position = com
        self.idx = node_idx
        self.root = Node(node_ndx, centre=com, size = 4*np.pi, leaf_num= 8)
        self.root.directs = np.array([[0., 0., 1.], #North
                                      [1., 0., 0.],[0., 1., 0.],[-1., 0., 0.],[0.,-1.,0.], # equator
                                      [0., 0., -1.0] # South
                                      ])
        for i in range(6):
            self.root.grids.append(Bitree(node_ndx + 'N%d'%(i), centre = 0.0, 
                                          size=np.pi, self.position, self.root.directs[i]))
        self.nodes = {}
        self.grids = {}
        self.iterateGrid()
        self.iterateNode()
        _grid_idx = [ [0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1], 
                      [5, 1, 2], [5, 2, 3], [5, 3, 4], [5, 4, 1]
                    ]
        for i in range(8):
            #three grids id in vertex A, B , C
            A, B, C = _grid_idx[i]
            idx = self.idx + str(i)
            directs = np.array([self.root.directs[A], 
				self.root.directs[B], 
				self.root.directs[C]])
            centre = np.sum(directs, axis=0)
            centre /= np.linalg.norm(centre)
            area = self._sphere_triang_area(directs[A], directs[B], directs[C])
            child = Node(idx, centre, area, leaf_num=4)
            Agrid = self.grepGrid(directs[A])
            Bgrid = self.grepGrid(directs[B])
            Cgrid = self.grepGrid(directs[C])
            child.grids = [Agrid,Bgrid,Cgrid ]
            child.directs = directs
            key = self.findChild(self.root, centre)
            self.root.children[key] = child 
        self.root.isLeafNode = False
        self.iterateGrid()
        self.iterateNode()
                
    def grepGrid(self, vector):
        """check if a grid(bitree) exists by distance of two vector
        
        """
        for idx, grid in self.grids.items():
            delta = np.linalg.norm(item.direct - vector)
            if delta < 0.001:
                return grid
        retrun None

    def subdivideNode(self, parent):
        parent.isLeafNode = False
        grids = parent.grids
        directs = parent.directs
        directs = np.array([directs[0],directs[1],directs[2],
                            directs[0] + directs[1],
                            directs[1] + directs[2],
                            directs[2] + directs[0]
                            ])
        #       0
        #     3   5
        #   1   4   2
        for i in range(3,6):
            grid = self.grepGrid(directs[i])
            if not grid:
                grid = Bitree(parent.idx + 'N%d'%(i), centre = 0.0, 
			      size=np.pi, self.position, directs[i])
            grids.append(grid)
        
        _grid_idx = [ [0,3,5], [3,1,4], [5,4,2], [3,4,5]]
        for i in range(4):
            A,B,C = _grid_idx[i]
            idx = parent.idx + str(i)
            child_directs = np.array([directs[A], directs[B], directs[C]])
	    centre = np.sum(child_directs, axis=0)		
	    centre /= np.linalg.norm(centre)
            area = self._sphere_triang_area(directs[A], directs[B], directs[C])
	    child = Node(idx, centre, area, leaf_num=4)
	    child.directs = child_directs
	    child.grids = [grids[A],grids[B],grids[C]]
	    parent.children[i] = child
	self.iterateGrid()
	self.iterateNode()

    def fill(self, conf_idx, value):
        """fill conf after generation 
        """
	grid_idx = conf_idx.split('C')[0]
	if grid_idx in self.grids:
            self.grids[grid_idx].fill(conf_idx, value)
        else:
            self.addNode(grid_idx)
            if grid_idx in self.grids:
		self.nodes[grid_idx].fill(conf_idx, value)
            else:
                raise Exception("Con't fill conf %s\n"%(conf_idx))

    def addNode(self,node_idx):
        node_idx = node_idx.split('N')[0]
        pre_idx, idxs = node_idx.split('R')
        pre_idx += 'R' 
        for i in idxs:
            pre_idx +=   str(i)
            if pre_idx not in self.nodes:
                self.subdivideNode(self.nodes[pre_idx[:-1]])
                self.iterateNode()

    def interpolation(self, vector, angle):
        neighbors = self.fi ndNeighbors(self.root, vector)
        v1 = neighbors[0].interpolation(angle)
        v2 = neighbors[1].interpolation(angle)
        v3 = neighbors[2].interpolation(angle)
        w1 = self._sphere_triang_area(neighbors[1].direct,neighbors[2].direct,vector)
        w2 = self._sphere_triang_area(neighbors[0].direct,neighbors[2].direct,vector)
        w3 = self._sphere_triang_area(neighbors[1].direct,neighbors[0].direct,vector)
        value = (v1 * w1 + v2 * w2 + v3 * w3)/(w1 + w2 + w3)
        return value

    def findNeighbors(self, node, vector):
        _neighbor = None
        if node == None:
            return None  
        elif node.isLeafNode:
            _neighbor = (node.grids[0], node.grids[1], node.grids[2])
            return _neighbor
        else:
            child = self.findChild(node, vector)
            return findNeighbors(node.children[child],vector)
            
    def findChild(self, node, vector):
	if node is self.root:
	    key = ''
	    for i in range(3):
		if vector[i] >= 0:
		    key += '+'
		else:
		    key += '-'
	    return DIRLOOKUP[key]
	else:
	    # three middle point of triangle is the directs of 4 th children
	    mids = node.children[3].directs
	    solve = np.linalg.solve(mids.T, vector)
	    if solve[0] < 0: return 2
	    if solve[1] < 0: return 0
	    if solve[2] < 0: return 1
	    return 3

    def iterateGrid(self):
        for conf in self._iterateGrid_help(self.root):
            if conf.idx not in self.grids:
                self.grids[conf.idx] = conf
        
    def _iterateGrid_help(self, node):
        """iterate all conf, not unique

        """
        for conf in node.grids:
            if conf != None:yield conf
        for child in node.children:
            for c in self._iterateGrid_help(child):
                yield c 

    def iterateNode(self):
        for n in self._iterateNode_hel p(self.root):
            if n.idx not in self.nodes:
                self.nodes[n.idx] = n

    def _iterateNode_help(self, node):
        if node != None:
            yield node 
        for child in node.children:
            for n in self._iterateNode_help(child):
                yield  n
    
    def _vet2ang(x, y):
        """get the angle of 2 vector

        """
        lx = np.sqrt(np.dot(x,x))
        ly = np.sqrt(np.dot(y,y))
        cos_angle = np.dot(x,y)/(lx * ly)
        angle = np.arccos(cos_angle)
        return angle

    def _sphere_triang_area(OA,OB,OC, r = 1):
        """get area of spherical triangle from 3 vectors (O point to surface).

        """
        a = vet2ang(OB,OC)
        b = vet2ang(OA,OC)
        c = vet2ang(OA,OB)
        cosA = (np.cos(a) - np.cos(b)*np.cos(c))/(np.sin(b)*np.sin(c))
        cosB = (np.cos(b) - np.cos(a)*np.cos(c))/(np.sin(a)*np.sin(c))
        cosC = (np.cos(c) - np.cos(b)*np.cos(a))/(np.sin(b)*np.sin(a))
        E = np.arccos(cosA) + np.arccos(cosB) + np.arccos(cosC) - np.pi
        return (E * r**2)



class Octree:
    """QuadTree to index sphere

    """
    def __init__(self, centre=(0.,0.,0.), size=12.0, symmetry=None):
        """init octree by symmetry, and add 
        
        com is centre of mass.
        node_idx will be like "T123R123N123C123"
        T:translocation, R:rotation, R:normal, C:configuration
        """
        self.nodes = {}
        self.grids = {}
        self.sym = symmetry #used by self.subdivideNode()
        self.root = Node('T', centre=centre, size = 12.0, leaf_num= 8)
        self.subdivideNode(self.root)
        if self.sym == 1:
            for key,i in DIRLOOKUP.items():
                if key[2]='-':
                    self.root.children[i]=None
        elif self.sym == 2:
            for key,i in DIRLOOKUP.items():
                if key[2]=='-' or key[1]=='-':
                    self.root.children[i]=None
        elif self.sym == 3:
            for key,i in DIRLOOKUP.items():
                if key[2]=='-' or key[1]=='-' or key[0]=='-':
                    self.root.children[i]=None
        else:
            if self.sym != None:
                raise Exception("translation symmetry has only 1,2 and 3.")
        #regeneration nodes and grids list for deleting Node.
        self.nodes = {} 
        self.grids = {}
        self.iterateGrid()
        self.iterateNode()
        
    def grepGrid(self, vector):
        """check if a grid(bitree) exists by distance of two vector
        
        """
        for idx, grid in self.grids.items():
            delta = np.linalg.norm(item.direct - vector)
            if delta < 0.001:
                return grid
        retrun None

    def _grid_positions(self,node,offset=None):
        if offset=None:offset = node.size
        offsets = np.array([[-offset,-offset,-offset],
                            [-offset,-offset,+offset],
                            [-offset,+offset,-offset],
                            [-offset,+offset,+offset],
                            [+offset,-offset,-offset],
                            [+offset,-offset,+offset],
                            [+offset,+offset,-offset],
                            [+offset,+offset,+offset]
                            ])
        return node.centre + offsets

    def _newCentre(self,node):
        offset = node.size/2.0
         return self._grid_positions(node,offset)
        
    def subdivideNode(self, parent):
        parent.isLeafNode = False
        newCentre = self._newCentre(parent)
        grid_positions = self._grid_positions(parent)
        for i in range(8):
             child = Node(parent.idx + str(i), centre=newCentre[i], 
                         size = parent.size/2.0, leaf_num= 8)
            for j, pos in enumerate(self._grid_positions(child)):
                 grid = self.grepGrid(pos)
                if not grid:
                    grid = Quadtree(pos, child.idx+'R%d'%(j))
            parent.children[i] = child
	self.iterateGrid()
	self.iterateNode()

    def fill(self, conf_idx, value):
        """fill conf after generation 
        """
	grid_idx = conf_idx.split('N')[0]
	if grid_idx in self.grids:
            self.grids[grid_idx].fill(conf_idx, value)
        else:
            self.addNode(grid_idx)
            if grid_idx in self.grids:
		self.grids[grid_idx].fill(conf_idx, value)
            else:
                raise  Exception("Con't fill conf %s\n"%(conf_idx))

    def addNode(self,node_idx):
        node_idx = node_idx.split('R')[0]
        pre_idx, idxs = node_idx.split('T')
        pre_idx += 'T' 
        for i in idxs:
            pre_idx += str(i)
            if pre_idx not in self.nodes:
                self.subdivideNode(self.nodes[pre_idx[:-1]])
                self.i terateNode()

    def interpolation(self, position, vector, angle):
        neighbors = self.findNeighbors(self.root, position)
        ndim = len(position)        
        v = np.zeros((8,ndim + 7))
        for i in range(8):
            v[i,0:ndim] = neighbors[i].centre
            v[i,ndim:] = neighbors[i].interpolation(vector, angle)
        for dim in range(ndim):
            vtx_delta = 2**(ndim - dim - 1)
            for vtx in range(vtx_delta):
                v[vtx, ndim:] += (  (v[vtx + vtx_delta, ndim:] - v[vtx, ndim:]) * 
                                (position[dim] - v[vtx,dim])/ (v[vtx + vtx_delta, dim] - v[vtx, dim])
                             )
        return v[0,ndim:]

    def findNeighbors(self, node, vector):
        _neighbor = None
        if node == None:
            return None  
        elif node.isLeafNode:
            _neighbor = node.grids
            return _neighbor
        else:
            child = self.findChild(node, vector)
            return findNeighbors(node.children[child],vector)
            
    def findChild(self, node, vector):
	key = ''
	for i in range(3):
	    if vector[i] >= node.centre[i]:
		key += '+'
	    else:
		key += '-'
	return DIRLOOKUP[key]

    def iterateGrid(self):
        for conf in self._iterateGrid_help(self.root):
            if conf.idx not in self.grids:
                self.g rids[conf.idx] = conf
        
    def _iterateGrid_help(self, node):
        """iterate all conf, not unique

        """
        for conf in node.grids:
            if conf != None:yield conf
        for child in node.children:
            for c in self._iterateGrid_help(child):
                yield c  

    def iterateNode(self):
        for n in self._iterateNode_hel p(self.root):
            self.nodes[n .idx] = n

    def _iterateNode_help(self, node):
        if node != None:
            yield node 
        for child in node.children:
            for n in self._iterateNode_help(child):
                yield  n 

    def _vet2ang(x, y):
        """get the  angle of 2 vector

        """
        lx = np.sqrt(np.dot(x,x))
        ly = np.sqrt(np.dot(y,y))
        cos_angle = np.dot(x,y)/(lx * ly)
        angle = np.arccos(cos_angle)
        return angle 

## ---------------------------------------------------------------------------------------------------##
class Grid:
    def __init__(self,centre=(0.,0.,0.), size=12.0, symmetry=None):
        self.mesh = Octree(centre=(0.,0.,0.), size=12.0, symmetry=None)
    
    def fill(conf_idx, value):
        self.mesh.fill(conf_idx, value)
    
    def interpolate(self, position, vector, angle):
        return self.mesh.interpolation(position, vector, angle)

    def load(self, filename):
        with open(filename,'r') as f:
            for i in f.readlines():
                if i[0] == '#':continue
                i = i.split()
                if len(i) != 7:continue
                idx = i[0]
                value = np.array([float(v) for v in i[1:]])
                self.mesh.fill(idx, value)
    
    def write(self, filename):
        with open(filename, 'w') as f:
            for conf in self._iter_conf(self):
                conf_str='%s'%(conf.idx) + ' %f'*7%tuple(conf.value) + '\n'
                f.write(conf_str)
    
    def gen_x(self):
        for conf in self._iter_conf(self):
            yield (conf.position, conf.vector, conf.angle)
    
    def _iter_conf(self):
        for quadtree in self.mesh.grids.values():
            for bitree in quadtree.grids.values():
                for conf in bitree.grids.values():
                    yield conf
    
            
        
if __name__ == "__main__":

    ### Object Insertion Test ###
    
    # So lets test the adding:
    import random
    import time

    #Dummy object class to test with
    class TestObject:
        def __init__(self, name, position):
            self.name = name
            self.position = position

    # Create a new octree, size of world
    myTree = Octree(40.0000)

    # Number of objects we intend to add.
    NUM_TEST_OBJECTS = 2000

    # Number of collisions we're going to test
    NUM_COLLISION_LOOKUPS = 2000

    # Insert some random objects and time it
    Start = time.time()
    for x in range(NUM_TEST_OBJECTS):
        name = "Node__" + str(x)
        pos = (random.uniform(-12.000, 12.000), random.uniform(-12, 12.00), random.uniform(-12.00, 12.00))
        testOb = TestObject(name, pos)
        myTree.insertNode(myTree.root, 12.000, myTree.root, testOb)
    End = time.time() - Start

    # print some results.
    print str(NUM_TEST_OBJECTS) + "-Node Tree Generated in " + str(End) + " Seconds"
    print "Tree Leaves contain a maximum of " + str(MAX_OBJECTS_PER_CUBE) + " objects each."

    ### Lookup Tests ###

    # Look up some random positions and time it
    Start = time.time()
    for x in range(NUM_COLLISION_LOOKUPS):
        pos = (random.uniform(-12.000, 12.000), random.uniform(-12.00, 12.00), random.uniform(-12.00, 12.00))
        result = myTree.findPosition(myTree.root, pos)
        
        ##################################################################################
        # This proves that results are being returned - but may result in a large printout
        # I'd just comment it out and trust me :)
        # print "Results for test at: " + str(pos)
        # if result != None:
        #    for i in result:
        #        print i.name, i.position,
        # print
        ##################################################################################
        
    End = time.time() - Start

    # print some results.
    print str(NUM_COLLISION_LOOKUPS) + " Collision Lookups performed in " + str(End) + " Seconds"
    print "Tree Leaves contain a maximum of " + str(MAX_OBJECTS_PER_CUBE) + " objects each."

    x = raw_input("Press any key (Wheres the any key?):")
