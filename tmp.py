import numpy as np
"""A Hierarchical Adaptive Meshmesh indexing tranlation and rotation space.

Combining octree (for indexing tranlation space) and quad-tree (for indexing sphere) to indexing all configurations of two rigid body.
"""
MAX_OBJECTS_PER_CUBE = 10


# This dictionary is used by the findBranch function, to return the correct branch index
DIRLOOKUP = {"3":0, "2":1, "-2":2, "-1":3, "1":4, "0":5, "-4":6, "-3":7}

#### End Globals ####

class conf:
    def __init__(self, node_idx, weight, value):
        self.idx = node_idx
        self.weight = weight
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
    def __init__(self,node_idx, centre, size=np.pi * 2, direction):
        self.direct = direction
        self.root = Node(node_idx, centre, size, leaf_num=2)
        self.root.grids.append(conf(node_idx+'C1',self.centre - size/2.0, 0.0))
        self.root.grids.append(conf(node_idx+'C0',self.centre, 0.0))
        self.root.grids.append(conf(node_idx+'C1',self.centre + size/2, 0.0))
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
        left.grids.append(conf(left.idx + 'C0', left.centre, 0.0))
        left.grids.append(parent.grids[1])
        right = parent.children[1]
        right.grids.append(parent.grids[1])
        right.grids.append(conf(left.idx + 'C0', right.centre, 0.0))
        right.grids.append(parent.grids[2])
    
    def fill(self, conf_idx, value):
        """fill conf after generation 
        """
        if conf_idx in self.grids:
            self.grids[conf_idx].value = value
        else:
            addNode(conf_idx)
            if conf_idx in self.grids:
                self.grids[conf_idx].value = value
            else:
                raise Exception("Con't fill conf %s\n"%(conf_idx))

    def addNode(self,node_idx):
        node_idx = node_idx.split('C')[0]
        pre_idx, idxs = node_idx.split('N')
        pre_idx += 'N' 
        for i in idxs:
            pre_idx +=  + str(i)
            if pre_idx not in self.nodes:
                self.subdivideNode(self.nodes[pre_idx[:-1]])
                self.iterateNode()

    def interpolation(self, angle):
        neighbors = self.findNeighbors(self.root, angle)
        v1 = neighbors[0].value
        v2 = neighbors[1].value
        w1 = angle - neighbors[0].weight
        w2 = neighbors[1].weight -angle
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
        for conf in node.grids:yield conf
        for child in node.children:
            for c in self._iterateGrid_help(child):
                yield c

    def iterateNode(self):
        for n in self._iterateNode_help(self.root):
            self.nodes[n.idx] = n

    def _iterateNode_help(self, node):
        yield node
        for child in node.children:
            yield child
            
class Quadtree:
    """QuadTree to index sphere

    """
    def __init__(self, com, node_idx):
        """init 6 grids and 8 triangles in root node.
        
        com is centre of mass.
        grid order is like '4', triangles child is anti o' clock from (1,1,1)
        """
        self.com = com
        self.idx = node_idx
        self.root = Node(node_ndx, centre=com, size = 4*np.pi, leaf_num= 8)
        _directs = np.array([[0., 0., 1.], #North
                            [1., 0., 0.],[0., 1., 0.],[-1., 0., 0.],[0.,-1.,0.], # equator
                            [0., 0., -1.0] # South
                            ])
        for i in range(6):
            self.root.grids.append(Bitree(node_ndx + 'N%d'%(i), centre = 0.0, size=np.pi, directs[i]))
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
            centre = _directs[A] + _directs[B] + _directs[C]
            centre /= np.linalg.norm(centre)
            area = self._sphere_triang_area(_directs[A], _directs[B], _directs[C])
            child = Node(idx, centre, area, leaf_num=4)
            Agrid = self.grepGrid(_directs[A])
            Bgrid = self.grepGrid(_directs[B])
            Cgrid = self.grepGrid(_directs[C])
            child.grids = [Agrid,Bgrid,Cgrid ]
            self.root.children[i] = child 
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
        diret = []
        for i in parent.grids:
            
        _centre = (parent.grids[])
        for i in ddrange(4):
            parent.children.append(Node(parent.idx+str(i),_centre[i], self.size/2.0, 2))
        left = parent.children[0]
        left.grids.append(parent.grids[0])
        left.grids.append(conf(left.idx + 'C0', left.centre, 0.0))
        left.grids.append(parent.grids[1])
        right = parent.children[1]
        right.grids.append(parent.grids[1])
        right.grids.append(conf(left.idx + 'C0', right.centre, 0.0))
        right.grids.append(parent.grids[2])
    
    def fill(self, conf_idx, value):
        """fill conf after generation 
        """
        if conf_idx in self.grids:
            grids[conf_idx].value = value
        else:
            addNode(conf_idx)
            if conf_idx in self.grids:
                self.grids[conf_idx].value = value
            else:
                raise Exception("Con't fill conf %s\n"%(conf_idx))

    def addNode(self,node_idx):
        node_idx = node_idx.split('C')[0]
        pre_idx, idxs = node_idx.split('N')
        pre_idx += 'N' 
        for i in idxs:
            pre_idx +=  + str(i)
            if pre_idx not in self.nodes:
                self.subdivideNode(self.nodes[pre_idx[:-1]])
                self.iterateNode()

    def interpolation(self, angle):
        neighbors = self.findNeighbors(self.root, angle)
        v1 = neighbors[0].value
        v2 = neighbors[1].value
        w1 = angle - neighbors[0].weight
        w2 = neighbors[1].weight -angle
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
        for conf in node.grids:yield conf
        for child in node.children:
            for c in self._iterateGrid_help(child):
                yield c

    def iterateNode(self):
        for n in self._iterateNode_help(self.root):
            self.nodes[n.idx] = n

    def _iterateNode_help(self, node):
        yield node
        for child in node.children:
            yield child
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





class OctNode:
    """OctNode store a cube space and its subspace.

    Args:
        centre(list): centre of cube
        size(float): side length of cube
        data(:obj:`list` | :obj:`tuple` | :obj:`str`, optional): some thing you want to save.
    
    Attributes:
        neighbors(list): the grids in the vertexs of cube. there is quad-tree on every grid.
        grids(list): the new gird in this node, for uniq loop.
        branches(list): the childs of this node.
    """
    def __init__(self, position, size, data=None):
        """init OctNode.

        """
        self.centre = centre
        self.size = size

        # All OctNodes will be leaf nodes at first
        # Then subdivided later as more objects get added
        self.isLeafNode = True

        # store our object, typically this will be one, but maybe more
        self.data = data
        
        self.branches = [None, None, None, None, None, None, None, None]
        
    def grids(self):
        x,y,z = self.centre
        s = self.size
        self._location = [[(x+s,y+s,z+s),None,None,None,None],
                          [None,None,None,None,None],
                          [None,None,None,None,None]]
        def neighors
        self.neighors = self.grids[ [self.grids[]],
                                    [] ]
        # The cube's bounding coordinates -- Not currently used
        self.ldb = (position[0] - (size / 2), position[1] - (size / 2), position[2] - (size / 2))
        self.ruf = (position[0] + (size / 2), position[1] + (size / 2), position[2] + (size / 2))
        


class Octree:
    def __init__(self, worldSize):
        # Init the world bounding root cube
        # all world geometry is inside this
        # it will first be created as a leaf node (ie, without branches)
        # this is because it has no objects, which is less than MAX_OBJECTS_PER_CUBE
        # if we insert more objects into it than MAX_OBJECTS_PER_CUBE, then it will subdivide itself.
        self.root = self.addNode((0,0,0), worldSize, [])
        self.worldSize = worldSize

    def addNode(self, position, size, objects):
        # This creates the actual OctNode itself.
        return OctNode(position, size, objects)

    def insertNode(self, root, size, parent, objData):
        if root == None:
            # we're inserting a single object, so if we reach an empty node, insert it here
            # Our new node will be a leaf with one object, our object
            # More may be added later, or the node maybe subdivided if too many are added
            # Find the Real Geometric centre point of our new node:
            # Found from the position of the parent node supplied in the arguments
            pos = parent.position
            # offset is halfway across the size allocated for this node
            offset = size / 2
            # find out which direction we're heading in
            branch = self.findBranch(parent, objData.position)
            # new center = parent position + (branch direction * offset)
            newCenter = (0,0,0)
            if branch == 0:
                # left down back
                newCenter = (pos[0] - offset, pos[1] - offset, pos[2] - offset )
                
            elif branch == 1:
                # left down forwards
                newCenter = (pos[0] - offset, pos[1] - offset, pos[2] + offset )
                
            elif branch == 2:
                # right down forwards
                newCenter = (pos[0] + offset, pos[1] - offset, pos[2] + offset )
                
            elif branch == 3:
                # right down back
                newCenter = (pos[0] + offset, pos[1] - offset, pos[2] - offset )

            elif branch == 4:
                # left up back
                newCenter = (pos[0] - offset, pos[1] + offset, pos[2] - offset )

            elif branch == 5:
                # left up forward
                newCenter = (pos[0] - offset, pos[1] + offset, pos[2] + offset )
                
            elif branch == 6:
                # right up forward
                newCenter = (pos[0] + offset, pos[1] - offset, pos[2] - offset )

            elif branch == 7:
                # right up back
                newCenter = (pos[0] + offset, pos[1] + offset, pos[2] - offset )
            # Now we know the centre point of the new node
            # we already know the size as supplied by the parent node
            # So create a new node at this position in the tree
            # print "Adding Node of size: " + str(size / 2) + " at " + str(newCenter)
            return self.addNode(newCenter, size, [objData])
        
        #else: are we not at our position, but not at a leaf node either
        elif root.position != objData.position and root.isLeafNode == False:
            
            # we're in an octNode still, we need to traverse further
            branch = self.findBranch(root, objData.position)
            # Find the new scale we working with
            newSize = root.size / 2
            # Perform the same operation on the appropriate branch recursively
            root.branches[branch] = self.insertNode(root.branches[branch], newSize, root, objData)
        # else, is this node a leaf node with objects already in it?
        elif root.isLeafNode:
            # We've reached a leaf node. This has no branches yet, but does hold
            # some objects, at the moment, this has to be less objects than MAX_OBJECTS_PER_CUBE
            # otherwise this would not be a leafNode (elementary my dear watson).
            # if we add the node to this branch will we be over the limit?
            if len(root.data) < MAX_OBJECTS_PER_CUBE:
                # No? then Add to the Node's list of objects and we're done
                root.data.append(objData)
                #return root
            elif len(root.data) == MAX_OBJECTS_PER_CUBE:
                # Adding this object to this leaf takes us over the limit
                # So we have to subdivide the leaf and redistribute the objects
                # on the new children. 
                # Add the new object to pre-existing list
                root.data.append(objData)
                # copy the list
                objList = root.data
                # Clear this node's data
                root.data = None
                # Its not a leaf node anymore
                root.isLeafNode = False
                # Calculate the size of the new children
                newSize = root.size / 2
                # distribute the objects on the new tree
                # print "Subdividing Node sized at: " + str(root.size) + " at " + str(root.position)
                for ob in objList:
                    branch = self.findBranch(root, ob.position)
                    root.branches[branch] = self.insertNode(root.branches[branch], newSize, root, ob)
        return root

    def findPosition(self, root, position):
        # Basic collision lookup that finds the leaf node containing the specified position
        # Returns the child objects of the leaf, or None if the leaf is empty or none
        if root == None:
            return None
        elif root.isLeafNode:
            return root.data
        else:
            branch = self.findBranch(root, position)
            return self.findPosition(root.branches[branch], position)
            

    def findBranch(self, root, position):
        # helper function
        # returns an index corresponding to a branch
        # pointing in the direction we want to go
        vec1 = root.position
        vec2 = position
        result = 0
        # Equation created by adding nodes with known branch directions
        # into the tree, and comparing results.
        # See DIRLOOKUP above for the corresponding return values and branch indices
        for i in range(3):
            if vec1[i] <= vec2[i]:
                result += (-4 / (i + 1) / 2)
            else:
                result += (4 / (i + 1) / 2)
        result = DIRLOOKUP[str(result)]
        return result

## ---------------------------------------------------------------------------------------------------##


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
