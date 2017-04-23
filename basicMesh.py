"""A Hierarchical Adaptive Mesh indexing tranlation and rotation space.

Combining octree (for indexing tranlation space) and 
cubed-sphere quaterion to index all configurations of two rigid body.

Octree grid order:like N\\N, this is a Z-like order
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
import numpy as np
import pdb
np.seterr(invalid='warn')

ORDER = 1
# This dictionary is used by the findBranch function, to return the correct branch index
DIRLOOKUP = {'+++':7, '-++':3, '--+':1, '+-+':5, '++-':6, '-+-':2, '---':0, '+--':4}
DIRLOOKUP4 = {'--':0, '-+':1, '+-':2, '++':3}
#### End Globals ####
#E_HIGH = 100.0
E_HIGH = 0.0 # default value set as 0.0 maybe more retionable
DIST_CUTOFF = 2.5 ** 2
HIGH = np.array([E_HIGH,E_HIGH,E_HIGH,E_HIGH,E_HIGH,E_HIGH,E_HIGH])
CONFS = [] # save confs in the order of creation
GOLDEN = (1 + 5 ** 0.5) / 2

class conf:
    def __init__(self, idx, location=None, q=None,values=HIGH):
        self.idx = idx
        self.pos =  {} # position in angle cubic
        self.loc = location # location  in  translation space (r, theta, phi)
        self.q = q
        self.values = values

class Node:
    """Node.

    """
    def __init__(self, idx=None, pos=None, size=None, leaf_num=None):
        self.error = 100.0
        self.pos = pos
        self.size = size
        self.isLeafNode  = True
        self.idx = idx
        self.leaf_num = leaf_num
        # children store region
        self.children = [None for i in range(self.leaf_num)] # the branches should have order
        # grids store mesh grids
        self.grids = []
        self.testgrid = []
        self.parent = None


class Octree:
    """QuadTree to index sphere

    """
    def __init__(self, idx, pos, size, next_level_id):
        """init octree by symmetry, and add 
        
        com is pos of mass.
        node_idx will be like "T123R123N123C123"
        T:translocation, R:rotation, R:normal, C:configuration
        """
        self.idx = idx
        self.pos = pos
        self.size = size
        self.nID = next_level_id
        self.allnodes = {}
        self.allgrids = {}
        self.gDict = dict()
        self.root = Node(idx, pos, size, leaf_num= 8)
        self.allnodes[self.root.idx] = self.root
        #print("init a Tree idx:%10s loc:r,%6.3f phi,%6.3f theta,%6.3f"%(self.idx, self.pos[0], self.pos[1],self.pos[2])  )
 
    def _triLinearInterp(self, pos, matrix):
        ndim = len(pos)        
        v = martix # position + values
        for dim in range(ndim):
            vtx_delta = 2**(ndim - dim - 1)
            for vtx in range(vtx_delta):
                v[vtx, ndim:] += (  (v[vtx + vtx_delta, ndim:] - v[vtx, ndim:]) * 
                                (pos[ dim] - v[vtx,dim])/ (v [vtx + vtx_delta, dim] - v[vtx, dim])
                               ) 
        return v[0,ndim:]

    def addNode(self, idx):
        node_idx = idx.split(self.nID)[0]
        if node_idx in self.allnodes:return
        #pdb.set_trace()
        #self._addNode_help(node_idx[:-1])
        idxs = ''
        while node_idx not in self.allnodes:
            idxs = node_idx[-1]+ idxs
            node_idx = node_idx[:-1]
        for i in idxs:
            if self.allnodes[node_idx].isLeafNode:
                self.subdivideNode(self.allnodes[node_idx])
            node_idx += i
   
    def _addNode_help(self, node_idx):
        #print(node_idx)
        if node_idx in self.allnodes:
            if self.allnodes[node_idx].isLeafNode:
                self.subdivideNode(self.allnodes[node_idx])
        else:
            self._addNode_help(node_idx[:-1])
            self.subdivideNode(self.allnodes[node_idx])

    def grepGrid(self, vector):
        """check if a grid(bitree) exists by distance of two vector
        print(self.root.children)
        
        """
        for grid in self.allgrids.values():
            delta = np.linalg.norm(grid.pos - vector)
            if delta < 0.001:
                return grid
        return None

    def _grid_positions(self,node,offset=None):
        if offset is None:offset = node.size
        offsets = np.array([[-offset[0],-offset[1],-offset[2]],
                            [-offset[0],-offset[1],+offset[2]],
                            [-offset[0],+offset[1],-offset[2]],
                            [-offset[0],+offset[1],+offset[2]],
                            [+offset[0],-offset[1],-offset[2]],
                            [+offset[0],-offset[1],+offset[2]],
                            [+offset[0],+offset[1],-offset[2]],
                            [+offset[0],+offset[1],+offset[2]]
                             ]) 
        return node.pos + offsets

    def _newCentre(self,node):
        offset = node.size/2.0
        return self._grid_positions(node,offset)
        
    def findNeighbors(self, node, vector):
        leaf = self.inWhichNode(node, vector)
        return leaf.grids

    def inWhichNode(self, node, vector):
        _neighbor = None
        if node == None:
            return None  
        elif node.isLeafNode:
            return node
        else:
            if len(vector) == 4 and node is self.root:
                child = np.abs(vector).argmax()
                vector /= vector[chlild]
                vector = np.arctan(vector)
                vector = np.delete(vector,child)
            else:
                child = self.findChild(node, vector)
            return self.inWhichNode(node.children[child],vector)
            
    def findChild(self, node, vector):
        key = ''
        for i in range(3):
            if vector[i] >= node.pos[i]:
                key += '+'
            else:
                key += '-'
        return DIRLOOKUP[key]


    def _np2str(self, q):
        npstr = ''
        for i in q:
            if i  < 0.0000001 and i> -0.0000001:i += 0.0000001
            npstr += 'N%.5f'%(i)
        return npstr


#########################################################################################
class Qtree(Octree):
    def __init__(self, idx, xyz):
        # {} is to save sphere ndx
        Octree.__init__(self, idx, {}, 16*3.14, 'C') #There 4 sphere in  root ^_^
        self.xyz = xyz
        self.leafNodes = set()
        self.leafNodes.add(self.root)
        self.subdivideNode(self.root)
        for i in range(4): self.subdivideNode(self.root.children[i])
        #self.fill(self.root.children[0].idx+'111111C7',0.0)

    def fill(self, idx, values):
        """fill conf after generation 
        """
        # import pdb 
        #pdb.set_trace()
        self.addNode(idx)
        # if idx:wtrR0Q2C7, grid_idx::wtrR0 or wtrR0Q2
        if idx in self.allgrids:
            self.allgrids[idx].values = values
        else:
            raise Exception("Con't fill conf %s\n"%(idx))
            
    def subdivideNode(self, parent):
        parent.isLeafNode = False
        self.leafNodes.remove(parent)
        childnum = 8
        if parent is self.root:
            childnum = 4
            for i in range(4):
                child = Node(self.idx + '%d'%(i), np.array([0.,0.,0.]), 
                         np.array([np.pi/4., np.pi/4., np.pi/4.]), 8)
                parent.children[i] = child
        else:
            newCentre = self._newCentre(parent)
            for i in range(8):
                child = Node(parent.idx + str(i), pos=newCentre[i],
                        size = parent.size/2.0, leaf_num= 8)
                parent.children[i] = child
        
        for i in range(childnum):
            child = parent.children[i]
            child.parent = parent
            self.leafNodes.add(child)
            self.allnodes[child.idx] = child
            for j, pos_ndx in enumerate(self._grid_positions(child)):
                idx = child.idx+'%s%d'%(self.nID,j)
                q = self._ndx2q(child.idx, pos_ndx)
                npstr = self._np2str(q)
                grid = None
                if npstr in self.gDict:
                    grid = self.gDict[npstr]
                else:
                    grid = conf(idx, self.xyz, q) 
                    CONFS.append(grid)
                    self.gDict[npstr] = grid
                    self.allgrids[grid.idx]=grid
                    # self.pos is translation space coord
                grid.pos[idx] = pos_ndx
                child.grids.append(grid)
            # add testgrid
        if parent is self.root: return
        idx = parent.idx+'%s8'%(self.nID)
        q = self._ndx2q(parent.idx, parent.pos)
        npstr = self._np2str(q)
        grid = None
        if npstr in self.gDict:
            grid = self.gDict[npstr]
        else:
            grid = conf(idx, self.xyz, q) 
            CONFS.append(grid)
            self.gDict[npstr] = grid
            self.allgrids[grid.idx]=grid
        parent.testgrid.append(grid)
        #self.iterateGrid()
        #self.iterateNode()
        #parent.testgrid.append(self.grepGrid(self._ndx2q(parent.idx, parent.pos)))
                
    def interpolation(self, q, node=None, neighbors=None):
        cubeID, ndx = self._q2ndx(q)
        if node is None: 
            node = self.root.children[cubeID]
        #else:
        #    pos = self.pos
        #    if np.dot(pos, pos) < DIST_CUTOFF
        #        return np.array([100.0,100.0,100.0,100.0,100.0,100.0,100.0])
        #    if np.dot(pos, pos) >  144.0:
        #        return np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0])
        if neighbors is None: 
            if not node.isLeafNode:
                node = self.inWhichNode(node, ndx)
            neighbors = node.grids
        if ORDER == 2 and  (None not in node.parent.children):
            return self._interpolation2(ndx, node)
        ndim = len(ndx)        
        v = np.zeros((8,ndim + 7))
        for i in range(8):
            v[i,0:ndim] = neighbors[i].pos[node.idx+'%s%d'%(self.nID,i)]# neigbor is conf, conf.pos is angle index
            v[i,ndim:] = neighbors[i].values
        for dim in range(ndim):
            vtx_delta = 2**(ndim - dim - 1)
            for vtx in range(vtx_delta):
                #if np.abs(v[vtx + vtx_delta, dim] - v[vtx, dim]) < 0.000001:
                   # pdb.set_trace()
                    #print(ndx[dim], v[vtx,dim],v[vtx + vtx_delta, dim])
                v[vtx, ndim:] += (  (v[vtx + vtx_delta, ndim:] - v[vtx, ndim:]) * 
                    (ndx[dim] - v[vtx,dim])/ (v[vtx + vtx_delta, dim] - v[vtx, dim])
                    )
        return v[0,ndim:]

    def _interpolation2(self, my_x,  node):
        children = node.parent.children
        grid_order = ((0,0),(0,1),(1,1),
                      (0,2),(0,3),(1,3),
                      (2,2),(2,3),(3,3),
                      
                      (0,4),(0,5),(1,5),
                      (0,6),(0,7),(1,7),
                      (2,6),(2,7),(3,7),
                      
                      (4,4),(4,5),(5,5),
                      (4,6),(4,7),(5,7),
                      (6,6),(6,7),(7,7)
                      )
        #   ___________
        #  1|1  3|1   3|3
        #   |  1 |  3  |
        #   |0__2|0___2|
        #   |1  3|1   3|
        #   |  0 |  2  |
        #  0|0__2|0___2|2
        #
        grids = []
        gridIdx = []
        for i, j in grid_order:
            grids.append(children[i].grids[j])
            gridIdx.append(children[i].idx + '%s%d'%(self.nID,j))
        grids = np.reshape(grids,(3,3,3))
        gridIdx = np.reshape(gridIdx,(3,3,3))
        x_dim1 = []
        y_dim1 = []
        for i in range(3):
            x_dim2 = []
            y_dim2 = []
            for j in range(3):
                x_dim3 = []
                y_dim3 = []
                for k in range(3):
                    x_dim3.append(grids[i][j][k].pos[gridIdx[i][j][k]][2])
                    y_dim3.append(grids[i][j][k].values)
                x_dim2.append(grids[i][j][0].pos[gridIdx[i][j][0]][1])
                y_dim2.append(self._interp_1D(x_dim3, y_dim3, my_x[2]))
            x_dim1.append(grids[i][0][0].pos[gridIdx[i][0][0]][0])
            y_dim1.append(self._interp_1D(x_dim2, y_dim2, my_x[1]))
        return self._interp_1D(x_dim1,y_dim1,my_x[0])
    
    def _interp_1D(self, xs, ys, my_x):
        if len(xs) == 1:
            return ys[0]
        if len(xs) == 2:
            x0, x1 = xs
            y0, y1 = ys
            return y0 + (my_x - x0) * (y1 - y0) / (x1 - x0)
        if len(xs) == 3:
            x0, x1, x2 = xs
            y0, y1, y2 = ys
            a1 = (y1 - y0) / (x1 - x0)
            a2 = (y2 - y0 - a1 * (x2 - x0)) / (x2 - x0) / (x2 - x1)
            return y0 + a1 * (my_x - x0) + a2 * (my_x - x0) * (my_x - x1)

    def _q2ndx(self, q):
        q = np.copy(q)
        cubeID = np.abs(q).argmax()
        q /= q[cubeID]
        q = np.arctan(q)
        ndx = np.delete(q, cubeID)
        return cubeID, ndx

    def _ndx2q(self, cubeID, ndx):
        cubeID = str(cubeID)
        cubeID = cubeID.replace(self.idx,'')[0]
        cubeID = int(cubeID)
        ndx = np.tan(ndx)
        q = np.insert(ndx, cubeID, 1.0)
        q /= np.linalg.norm(q)
        if q[0] < 0:
            q = -q
        if q[1] < 0: #water symmstery
            q[0], q[1], q[2], q[3] = -q[1], q[0], q[3], -q[2]
        if q[0]== q[1] == 0:
            if q[2] ==  0 or q[3] ==  0:
                q = np.array([0., 0., 1., 0.])
        if q[2]== q[3] == 0:
            if q[0] ==  0 or q[1] ==  0:
                q = np.array([1., 0., 0., 0.])
        return q


class Sphere:
    def __init__(self, idx, r, next_level_id='Q'):
        self.r = r
        self.idx = idx
        self.nID = next_level_id
        self.allnodes = {}
        self.allgrids = {}
        self.gDict = dict()
        self.root = Node(idx, r, 4*np.pi*r, leaf_num= 6)
        self.allnodes[self.root.idx] = self.root
        self.leafNodes = set()
        self.leafNodes.add(self.root)
        self.subdivideNode(self.root)
    
    def fill(self, idx, values):
        """fill conf after generation 
        """
        #import pdb 
        #pdb.set_trace()
        self.addNode(idx)
        # if idx:wtrR0Q2C7, grid_idx::wtrR0 or wtrR0Q2
        grid_idx = idx[ :idx.find(self.nID) + 2]
        if grid_idx in self.allgrids:
            self.allgrids[grid_idx].fill(idx, values)
        else:
            raise Exception("Con't fill conf %s\n"%(idx))

    def subdivideNode(self, parent):
        parent.isLeafNode = False
        self.leafNodes.remove(parent)
        if parent is self.root:
            # +x, +y, +z, -x, -y, -z 
            #signs = [1., 1., 1., -1., -1.,-1. ]
            num = 4 # sym = 2
            #num = 6 # no sym 
            for i in range(num):
                child = Node(parent.idx + str(i), np.array([0., 0.]), 
                             np.array([np.pi/4., np.pi/4.]), 4)
                parent.children[i] = child
                self.leafNodes.add(child)
                self.allnodes[child.idx] = child
                child.parent = parent
                self.subdivideNode(child)
            return
        else:
            newCentre = self._newCentre(parent)
            for i in range(4):
                xyz = self._ndx2xyz(parent.idx, newCentre[i])
                if xyz[1] <  -0.00001 or xyz[2] <  -0.00001: continue
                child = Node(parent.idx + str(i), pos=newCentre[i],
                             size = parent.size/2.0, leaf_num= 4)
                parent.children[i] = child
                self.leafNodes.add(child)
                self.allnodes[child.idx] = child

        for i in range(4):
            child = parent.children[i]
            if child is None: continue
            child.parent = parent
            for j, ndx in enumerate(self._grid_positions(child)):
                idx = child.idx+'%s%d'%(self.nID,j)
                xyz = self._ndx2xyz(child.idx, ndx)
                npstr = self._np2str(xyz)
                grid = None
                if npstr in self.gDict:
                    grid = self.gDict[npstr]
                else:
                    grid = Qtree(idx, xyz)
                    self.gDict[npstr] = grid
                    self.allgrids[grid.idx]=grid
                grid.pos[idx] = ndx
                child.grids.append(grid)
        idx = parent.idx+'%s8'%(self.nID)
        xyz = self._ndx2xyz(parent.idx, parent.pos)
        npstr = self._np2str(xyz)
        grid = None
        if npstr in self.gDict:
            grid = self.gDict[npstr]
        else:
            grid = Qtree(idx, xyz)
            self.gDict[npstr] = grid
            self.allgrids[grid.idx]=grid
        parent.testgrid.append(grid)
    
    def interpolation(self, xyz, q, node=None, neighbors=None):
        cubeID, ndx = self._xyz2ndx(xyz)
        if node is None:
            node = self.root.children[cubeID]
        if neighbors is None:
            if not node.isLeafNode:
                node = self.inWhichNode(node, ndx)
        if ORDER == 1:
            return self._interpolation1(ndx, q, node)
        if ORDER == 2:
            if None in node.parent.children:
                return self._interpolation1(ndx, q, node)
            return self._interpolation2(ndx, q, node)
    
    def _interpolation1(self, my_x, q, node):
        neighbors1 = (node.grids[:2], node.grids[2:])
        gridIdx  = [node.idx+'%s%d'%(self.nID,i) for i in range(4)]
        gridIdx = (gridIdx[:2], gridIdx[2:])
        x_dim1 = []
        y_dim1 = []
        for i in range(2):
            x_dim2 = []
            y_dim2 = []
            for j in range(2):
                x_dim2.append(neighbors1[i][j].pos[gridIdx[i][j]][1])
                y_dim2.append(neighbors1[i][j].interpolation(q))
            x_dim1.append(neighbors1[i][0].pos[gridIdx[i][0]][0])
            y_dim1.append(self._interp_1D(x_dim2, y_dim2, my_x[1]))
        return self._interp_1D(x_dim1,y_dim1,my_x[0])
                
    def _interpolation2(self, my_x, q, node):
        children = node.parent.children
        grid_order = ((0,0),(0,1),(1,1),
                      (0,2),(0,3),(1,3),
                      (2,2),(2,3),(3,3))
        #   ___________
        #  1|1  3|1   3|3
        #   |  1 |  3  |
        #   |0__2|0___2|
        #   |1  3|1   3|
        #   |  0 |  2  |
        #  0|0__2|0___2|2
        #
        grids = []
        gridIdx = []
        for i, j in grid_order:
            grids.append(children[i].grids[j])
            gridIdx.append(children[i].idx + '%s%d'%(self.nID,j))
        grids = (grids[:3],grids[3:6],grids[6:])
        gridIdx = (gridIdx[:3],gridIdx[3:6],gridIdx[6:])
        x_dim1 = []
        y_dim1 = []
        for i in range(3):
            x_dim2 = []
            y_dim2 = []
            for j in range(3):
                x_dim2.append(grids[i][j].pos[gridIdx[i][j]][1])
                y_dim2.append(grids[i][j].interpolation(q))
            x_dim1.append(grids[i][0].pos[gridIdx[i][0]][0])
            y_dim1.append(self._interp_1D(x_dim2, y_dim2, my_x[1]))
        return self._interp_1D(x_dim1,y_dim1,my_x[0])
    
    def _interp_1D(self, xs, ys, my_x):
        if len(xs) == 1:
            return ys[0]
        if len(xs) == 2:
            x0, x1 = xs
            y0, y1 = ys
            return y0 + (my_x - x0) * (y1 - y0) / (x1 - x0)
        if len(xs) == 3:
            x0, x1, x2 = xs
            y0, y1, y2 = ys
            a1 = (y1 - y0) / (x1 - x0)
            a2 = (y2 - y0 - a1 * (x2 - x0)) / (x2 - x0) / (x2 - x1)
            return y0 + a1 * (my_x - x0) + a2 * (my_x - x0) * (my_x - x1)
                
    def addNode(self, idx):
        node_idx = idx.split(self.nID)[0]
        if node_idx in self.allnodes:return
        idxs = ''
        while node_idx not in self.allnodes:
            idxs = node_idx[-1]+ idxs
            node_idx = node_idx[:-1]
        for i in idxs:
            if self.allnodes[node_idx].isLeafNode:
                self.subdivideNode(self.allnodes[node_idx])
            node_idx += i

    def _grid_positions(self,node,offset=None):
        if offset is None:offset = node.size

        offsets = np.array([[-offset[0],-offset[1]],
                            [-offset[0],+offset[1]],
                            [+offset[0],-offset[1]],
                            [+offset[0],+offset[1]]
                            ]) 
        return node.pos + offsets

    def _newCentre(self,node):
        offset = node.size/2.0
        return self._grid_positions(node,offset)
        
    def findNeighbors(self, node, vector):
        leaf = self.inWhichNode(node, vector)
        return leaf.grids

    def inWhichNode(self, node, vector):
        _neighbor = None
        if node == None:
            return None  
        elif node.isLeafNode:
            return node
        else:
            if len(vector) == 3 and node is self.root:
                child, vector = self._xyz2ndx
                node = self.root
            else:
                child = self.findChild(node, vector)
            return self.inWhichNode(node.children[child],vector)
            
    def findChild(self, node, vector):
        key = ''
        for i in range(2):
            if vector[i] >= node.pos[i]:
                key += '+'
            else:
                key += '-'
        return DIRLOOKUP4[key]

    def _ndx2xyz(self, cubeID, ndx):
        cubeID = str(cubeID)
        cubeID = cubeID.replace(self.idx,'')[0]
        cubeID = int(cubeID)
        ndx = np.tan(ndx)
        sign = 1.
        if cubeID > 2:sign = -sign
        cubeID %= 3
        _xyz = np.insert(ndx, cubeID, sign)
        _xyz /= np.linalg.norm(_xyz)
        _xyz *= self.r
        return _xyz
    
    def _xyz2ndx(self, xyz, r = None):
        _xyz = np.copy(xyz)
        _xyz /= np.linalg.norm(_xyz)
        cubeID = np.abs(_xyz).argmax()
        if _xyz[cubeID] < 0.:
            _xyz /=  - _xyz[cubeID]
            ndx = np.delete(_xyz, cubeID)
            cubeID += 3
        else:
            _xyz /=  _xyz[cubeID]
            ndx = np.delete(_xyz, cubeID)
        ndx = np.arctan(ndx)
        return cubeID, ndx

    def _np2str(self, q):
        npstr = ''
        for i in q:
            if i  < 0.0000001 and i> -0.0000001:i += 0.0000001
            npstr += 'N%.5f'%(i)
        return npstr

class Rtree:
    def __init__(self, idx, pos, size, next_level_id='S'):
        self.idx = idx
        self.pos = pos
        self.size = size
        self.nID = next_level_id
        self.allnodes = {}
        self.allgrids = {}
        self.gDict = dict()
        self.root = Node(idx, pos, size, leaf_num= 2)
        self.allnodes[self.root.idx] = self.root
        self.leafNodes = set()
        self.leafNodes.add(self.root)
        self.subdivideNode(self.root)
        for child in self.root.children:
            self.subdivideNode(child)
    
    def fill(self, idx, values):
        """fill conf after generation 
        """
        #import pdb 
        #pdb.set_trace()
        self.addNode(idx)
        # if idx:wtrR0Q2C7, grid_idx::wtrR0 or wtrR0Q2
        grid_idx = idx[ :idx.find(self.nID) + 2]
        if grid_idx in self.allgrids:
            self.allgrids[grid_idx].fill(idx, values)
        else:
            raise Exception("Con't fill conf %s\n"%(idx))

    def subdivideNode(self, parent):
        parent.isLeafNode = False
        self.leafNodes.remove(parent)
        newCentre = self._newCentre(parent)
        for i in range(2):
            child = Node(parent.idx + str(i), pos=newCentre[i], 
                         size = parent.size/2.0, leaf_num= 2)
            parent.children[i] = child
            child.parent = parent
            self.leafNodes.add(child)
            self.allnodes[child.idx] = child
            for j, r in enumerate(self._grid_positions(child)):
                idx = child.idx+'%s%d'%(self.nID,j)
                npstr = self._np2str(r)
                grid = None
                if npstr in self.gDict:
                    grid = self.gDict[npstr]
                else:
                    grid = Sphere(idx, r)
                    self.gDict[npstr] = grid
                    self.allgrids[grid.idx]=grid
                child.grids.append(grid)
        idx = parent.idx+'%s8'%(self.nID)
        r = parent.pos
        npstr = self._np2str(r)
        if npstr in self.gDict:
            grid = self.gDict[npstr]
        else:
            grid = Sphere(idx, r)
            self.gDict[npstr] = grid
            self.allgrids[grid.idx]=grid
        parent.testgrid.append(grid)

    def interpolation(self, xyz, q, node=None, neighbors=None):
        if node is None:
            node = self.root
        r = np.linalg.norm(xyz)
        if neighbors is None:
            if not node.isLeafNode:
                node = self.inWhichNode(node, r)
        if ORDER == 1:
            neighbors = node.grids
        if ORDER == 2:
            children = node.parent.children
            neighbors = (children[0].grids[0], children[0].grids[1], children[1].grids[1])
        xs = []
        ys = []
        for n in neighbors:
            xs.append(n.r)
            ys.append(n.interpolation(xyz, q))
        my_y = self._interp_1D(xs, ys, r)
        return my_y

    def _interp_1D(self, xs, ys, my_x):
        if len(xs) == 1:
            return ys[0]
        if len(xs) == 2:
            x0, x1 = xs
            y0, y1 = ys
            return y0 + (my_x - x0) * (y1 - y0) / (x1 - x0)
        if len(xs) == 3:
            x0, x1, x2 = xs
            y0, y1, y2 = ys
            a1 = (y1 - y0) / (x1 - x0)
            a2 = (y2 - y0 - a1 * (x2 - x0)) / (x2 - x0) / (x2 - x1)
            return y0 + a1 * (my_x - x0) + a2 * (my_x - x0) * (my_x - x1)


    def addNode(self, idx):
        node_idx = idx.split(self.nID)[0]
        if node_idx in self.allnodes:return
        idxs = ''
        while node_idx not in self.allnodes:
            idxs = node_idx[-1]+ idxs
            node_idx = node_idx[:-1]
        for i in idxs:
            if self.allnodes[node_idx].isLeafNode:
                self.subdivideNode(self.allnodes[node_idx])
            node_idx += i

    def _grid_positions(self,node,offset=None):
        if offset is None:offset = node.size
        offsets = np.array([-offset,+offset]) 
        return node.pos + offsets

    def _newCentre(self,node):
        offset = node.size/2.0
        return self._grid_positions(node,offset)
        
    def findNeighbors(self, node, vector):
        leaf = self.inWhichNode(node, vector)
        return leaf.grids

    def inWhichNode(self, node, vector):
        if node == None:
            return None  
        elif node.isLeafNode:
            return node
        else:
            child = self.findChild(node, vector)
            return self.inWhichNode(node.children[child],vector)
            
    def findChild(self, node, vector):
        if vector >= node.pos:
            return 1
        else:
            return 0

    def _np2str(self, i):
        npstr = ''
        if i  < 0.0000001 and i> -0.0000001:i += 0.0000001
        npstr += 'N%.5f'%(i)
        return npstr




################### Above: basic tree, node, grid(conf) class ######################
################### Below: basic mesh class ########################################
import pickle
class  mesh:
    def __init__(self):
        self.mesh = Rtree('wtr_wtrR', 7.0, 5.0)
        self.confs = set()
        self.confs.update(self._iter_conf())
        self.n = len(self.confs)
    def updateDatabase(self, filename):
        #self.confs = set()
        self.mesh = pickle.load(open(filename, "rb"))
        self.confs.update(self._iter_conf())
        self.n = len(self.confs)
        print("There are %d confs in this mesh."%(self.n))
        
    def dumpDatabase(self, filename):
        pickle.dump(self.mesh, open(filename, "wb"), True)
        print("Total %d confs saved in database %s"%(self.n, filename))

    def fill(self, conf_idx, values):
        self.mesh.fill(conf_idx, values)
    
    def fill_with_f(self, f, confs = None,filename='mesh.dat'):
        if confs is None: confs = self.confs
        n = len(confs)
        n_count = 1
        for conf in confs:
            if n_count%1000 == 1:
                print('-'*8 + "filling %10d/%d"%(n_count,n)+'-'*8)
            n_count += 1
            #self.fill(conf.idx, f(conf.loc, conf.q))
            conf.values =  f(conf.loc, conf.q)
        #self.save(confs)
        
    def interpolate(self, X, q):
        return self.mesh.interpolation(X, q)

    def load(self, filename):
        load_count = 0
        with open(filename,'r') as f:
            for i in f:
                if load_count%10000 == 0: print('>'*8+"Loaded %10d conf"%(load_count)+'<'*8)
                load_count += 1
                if i[0] == '#':continue
                i = i.split()
                if len(i) != 8:continue
                idx = i[0]
                #if idx=="wtr_wtrR360Q60C0":pdb.set_trace()
                values = np.array([float(v) for v in i[1:]])
                self.mesh.fill(idx, values)
        self._iter_conf()
    
    def save(self, confs=None, filename='mesh.dat'):
        if confs is None:confs=CONFS
        if self.database_name != None: filename = self.database_name
        with open(filename, 'w') as f:
            for conf in confs:
                conf_str='%s'%(conf.idx) + ' %f'*7%tuple(conf.values) + '\n'
                f.write(conf_str)

    def gen_x(self):
        for conf in self._iter_conf(self):
            yield (conf.loc, conf.q)
    
    def _iter_conf(self, node=None):
        if hasattr(node, 'q'): # type is conf
            yield node
        if node is None: node = self.mesh
        if hasattr(node, 'gDict'):
            for n in node.gDict.values():
                for c in self._iter_conf(n):
                    yield c
        if hasattr(node, '__len__'):
            for n in node:
                for c in self._iter_conf(n):
                    yield c

    def _grid_conf(self):
        for Qtree in self.Qtrees():
            yield Qtree.root.children[0].grids[0]

    def QLeafNodes(self):
        for Qtree in self.Qtrees():
            for n in Qtree.leafNodes:
                n.tree = Qtree
                yield n

    def SLeafNodes(self):
        for Stree in self.Strees():
            for n in Stree.leafNodes:
                n.tree = Stree
                yield n

    def RLeafNodes(self):
        for n in self.mesh.leafNodes:
            n.tree = self.mesh
            yield n

    def Qtrees(self):
        for Stree in self.mesh.gDict.values():
            for Qtree in Stree.gDict.values():
                yield Qtree
            
    def Strees(self):
        return self.mesh.gDict.values()
        
    
        
            
        
        
if __name__ == "__main__":
    pass
