import numpy as np
import pdb
np.seterr(invalid='warn')
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
#E_HIGH = 100.0
E_HIGH = 0.0 # default value set as 0.0 maybe more retionable
DIST_CUTOFF = 2.5 ** 2
HIGH = np.array([E_HIGH,E_HIGH,E_HIGH,E_HIGH,E_HIGH,E_HIGH,E_HIGH])

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
        self.root = Node(idx, pos, size, leaf_num= 8)
        #print("init a Tree idx:%10s loc:r,%6.3f phi,%6.3f theta,%6.3f"%(self.idx, self.pos[0], self.pos[1],self.pos[2])  )

#    def fill(self, idx, values):
#         """fill conf after generation 
#         """
#         #import pdb 
#         #pdb.set_trace()
#         self.addNode(idx)
#         # if idx:wtrR0Q2C7, grid_idx::wtrR0 or wtrR0Q2
#         grid_idx = idx[ :idx.find(self.nID) + 2]
#         if grid_idx in self.allgrids:
#            if hasattr(self.allgrids[grid_idx], 'fill'):
#                self.allgrids[grid_idx].fill(idx, values)
#            else:
#                self.allgrids[idx].values = values
#         else:
#             raise Exception("Con't fill conf %s\n"%(idx))
 
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
         grid_idx = idx[ :idx.find(self.nID) + 2]
         idxs = node_idx.replace(self.idx, '')
         pre_idx = self.idx
         for i in idxs:
             pre_idx  +=   str(i)
             if grid_idx not in self.allgrids:
                 self.subdivideNode(self.allnodes[pre_idx[:-1]])

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

    def iterateGrid(self):
        self.allgrids = {}
        for conf in self._iterateGrid_help(self.root):
            if conf.idx not in self.allgrids:
                self.allgrids[conf.idx] = conf
        
    def _iterateGrid_help(self, node):
        """iterate all conf, not unique

        """
        for conf in node.grids:
            if conf != None:
                yield conf
        for conf in node.testgrid:
            if conf != None:
                yield conf
                #print(conf.idx,conf)
        for child in node.children:
            if child is None:continue
            for c in self._iterateGrid_help(child):
                yield c  

    def iterateNode(self):
        self.allnodes = {}
        self.leafNodes = {}
        for n in self._iterateNode_help(self.root):
            if n.idx not in self.allnodes:
                self.allnodes[n.idx] = n
                if n.isLeafNode is True:
                    self.leafNodes[n.idx] = n 

    def _iterateNode_help(self, node):
        if node != None:
            yield node 
        for child in node.children:
            if child == None:continue
            for n in self._iterateNode_help(child):
                yield  n 


## ---------------------------------------------------------------------------------------------------##
class Qtree(Octree):
    def __init__(self, idx, pos):
        Octree.__init__(self, idx, pos, 16*3.14, 'C') #There 4 sphere in  root ^_^
        self.root.isLeafNode = False
        for i in range(4):
            #if "wtr_wtrR23Q12" == self.idx + str(i):pdb.set_trace()
            child = Node(self.idx + '%d'%(i), np.array([0.,0.,0.]), 
                         np.array([np.pi/4., np.pi/4., np.pi/4.]), 8)
            for j, pos_ndx in enumerate(self._grid_positions(child)):
                idx = child.idx+'%s%d'%(self.nID,j)
                q = self._ndx2q(child.idx, pos_ndx)
                grid = self.grepGrid(q)
                if not grid:
                    grid = conf(idx, self.pos, q)
                grid.pos[idx] = pos_ndx
                child.grids.append(grid)
            testidx = child.idx+'%d%s7'%(i, self.nID)
            q = self._ndx2q(child.idx, child.pos)
            grid = self.grepGrid(q)
            if not grid:
                grid = conf(testidx, self.pos, q)
            child.testgrid.append(grid)
            child.parent = self.root
            self.root.children[i] =  child
            self.iterateGrid()
            self.iterateNode()
            #self.subdivideNode(child)

    def fill(self, idx, values):
        """fill conf after generation 
        """
        #import pdb 
        #pdb.set_trace()
        self.addNode(idx)
        # if idx:wtrR0Q2C7, grid_idx::wtrR0 or wtrR0Q2
        if idx in self.allgrids:
            self.allgrids[idx].values = values
        else:
            raise Exception("Con't fill conf %s\n"%(idx))

#    def grepGrid(self, q):
#        for grid in self.allgrids.values():
#            delta1 = np.linalg.norm(grid.q - q)
#            delta2 = np.linalg.norm(grid.q + q)
#            if delta1 < 0.0001 or delta2 < 0.0001:
#                return grid
#            # For water only
#            q[0], q[1], q[2], q[3] = -q[1], q[0], q[3], -q[2]
#
#        return None
    
    def grepGrid(self, q):
        # For water like mol only
        p =np.zeros(4)
        p[0], p[1], p[2], p[3] = -q[1], q[0], q[3], -q[2]
        for grid in self.allgrids.values():
            delta1 = np.linalg.norm(grid.q - q)
            delta2 = np.linalg.norm(grid.q + q)
            delta3 = np.linalg.norm(grid.q + p)
            delta4 = np.linalg.norm(grid.q - p)
            if delta1 < 0.0001 or delta2 < 0.0001 or delta3 < 0.0001 or delta4 < 0.001:
                return grid
        return None

    def subdivideNode(self, parent):
        parent.isLeafNode = False
        newCentre = self._newCentre(parent)
        for i in range(8):
            child = Node(parent.idx + str(i), pos=newCentre[i], 
                         size = parent.size/2.0, leaf_num= 8)
            for j, pos_ndx in enumerate(self._grid_positions(child)):
                idx = child.idx+'%s%d'%(self.nID,j)
                q = self._ndx2q(child.idx, pos_ndx)
                grid = self.grepGrid(q)
                if not grid:
                    # self.pos is translation space coord
                    grid = conf(idx, self.pos, q) 
                grid.pos[idx] = pos_ndx
                child.grids.append(grid)
            testidx = child.idx+'%d%s7'%(i,self.nID)
            q = self._ndx2q(child.idx, child.pos)
            grid = self.grepGrid(q)
            if not grid:
                grid = conf(testidx, self.pos, q)
            child.testgrid.append(grid)
            parent.children[i] = child
            child.parent = parent
        self.iterateGrid()
        self.iterateNode()
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

    def _q2ndx(self, q):
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
        return q

class Rtree(Octree):
    def __init__(self, idx, symmetry=2):
        Octree.__init__(self,idx,
                        np.array([7., np.pi/4., np.pi/2.]), # idx,pos(r, phi(x), theta(xy))
                        np.array([5, np.pi/4., np.pi/2.]), 'Q') # size and rotation level idx
        self.subdivideNode(self.root)
        for i in self.root.children:
            self.subdivideNode(i)
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
    def interpolation(self, pos, q, node=None, neighbors=None):
        #if np.dot(pos, pos) < DIST_CUTOFF:
        #    return np.array([100.0,100.0,100.0,100.0,100.0,100.0,100.0])
        #if np.dot(pos, pos) >  144.0:
        #    return np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0])
        if node is None: node = self.root
        if neighbors is None: neighbors = self.findNeighbors(node, pos)
        ndim = len(pos)        
        v = np.zeros((8,ndim + 7))
        for i in range(8):
            v[i,0:ndim] = neighbors[i].pos
            v[i,ndim:] = neighbors[i].interpolation(q)
        for dim in range(ndim):
            vtx_delta = 2**(ndim - dim - 1)
            for vtx in range(vtx_delta):
                if v [vtx + vtx_delta, dim] ==  v[vtx, dim]:continue # avoid theta = pi/2, phi1=phi2
                v[vtx, ndim:] += ((v[vtx + vtx_delta, ndim:] - v[vtx, ndim:]) * 
                    (pos[dim] - v[vtx,dim])/ (v [vtx + vtx_delta, dim] - v[vtx, dim])
                    )
        return v[0,ndim:]
    
    def subdivideNode(self, parent):
        parent.isLeafNode = False
        newCentre = self._newCentre(parent)
        for i in range(8):
            child = Node(parent.idx + str(i), pos=newCentre[i], 
                         size = parent.size/2.0, leaf_num= 8)
            for j, pos in enumerate(self._grid_positions(child)):
                grid = self.grepGrid(pos)
                if not grid:
                    idx = child.idx+'%s%d'%(self.nID,j)
                    grid = Qtree(idx, pos)
                child.grids.append(grid)
            child.testgrid.append(Qtree(child.idx+'%d%s7'%(i, self.nID), child.pos))
            parent.children[i] = child
            child.parent = parent
        self.iterateGrid()
        self.iterateNode()
        parent.testgrid.append(self.grepGrid(parent.pos))

import copy
import pickle
class Grid:
    def __init__(self):
        self.mesh = Rtree('wtr_wtrR')
        self.confs = set()
        self.confs.update(self._iter_conf())
        self.n = len(self.confs)
    def updateDatabase(self, filename):
        self.mesh = pickle.load(open(filename, "rb"))
        self.confs.update(self._iter_conf())
        self.n = len(self.confs)
        print("There are %d confs in this mesh."%(self.n))
        
    def dumpDatabase(self, filename):
        pickle.dump(self.mesh, open(filename, "wb"), True)
        print("Total %d confs saved in database %s"%(self.n, filename))

    def fill(self, conf_idx, values):
        self.mesh.fill(conf_idx, values)
    def refine(self, f, toPDB, err_cutoff=1.0, filename='mesh.dat'):
        self.confs.update(self._iter_conf())
        self.n = len(self.confs)
        self.database_name = filename # set global filename for self.save
        self.Q_refine_count = 0
        self.R_refine_count = 0
        self.max_err_conf = None
        self.toPDB = toPDB
        self.f = f
        self.fill_with_f(self.f)
        self.E_CUTOFF = err_cutoff
        self._refine_Q()
        self._refine_R()
        self.n = len(self.confs)
        self.database_name = None
    
    def fill_with_f(self, f, confs = None,filename='mesh.dat'):
        if confs is None: confs = self.confs
        n = len(confs)
        n_count = 1
        for conf in confs:
            if n_count%1000 == 1:
                print('-'*8 + "filling %10d/%d"%(n_count,n)+'-'*8)
            n_count += 1
            self.fill(conf.idx, f(conf.loc, conf.q))
        #self.save(confs)
    

    def _refine_R(self):
        oldconfs = copy.copy(self.confs)
        fine = False
        while not fine:
            print("\n%d th time translocation R refinement"%(self.R_refine_count))
            self.R_refine_count += 1
            fine = True
            for leaf in self.RLeafNodes():
                if leaf.error < self.E_CUTOFF: continue
                # I want to restrict the density in very close
                # if leaf.pos[0] <  3 and leaf.size[0] < 0.5: continue
                node = leaf
                tree = leaf.tree
                node_err = 0
                testgrids = set()
                for Qtree in node.testgrid:
                    for conf in Qtree.allgrids.values():
                        testgrids.add(conf)
                for g in testgrids:
                    g_iterp = tree.interpolation(g.loc, g.q, node=node, neighbors=node.grids)
                    err = np.abs(g_iterp - g.values)
                    err = err[0]
                    if err > node_err:
                        node_err = err
                        self.max_err_conf = g
                if node_err < node.error: node.error = node_err
                if node.error > self.E_CUTOFF:
                    printStr=("max  error  is %5.2f\n"%(node.error)+
                                  "conf:%15s"%(self.max_err_conf.idx)+
                              " %5.2f"*3%tuple(self.max_err_conf.loc)+
                              " %5.2f"*4%tuple(self.max_err_conf.q) +
                              '\n' +
                              "size of node %6.3f\n"%(node.size[0]) +
                              "conf values is " +
                              " %5.2f"*7%tuple(self.max_err_conf.values)
                              )  
                    print(printStr)
                    tree.subdivideNode(node)
                    fine = False 
            self.confs.update(self._iter_conf())
            newconfs = self.confs.difference(oldconfs)
            print("total %15d confs, %15d new confs\n"%(len(self.confs),len(newconfs)))
            oldconfs = copy.copy(self.confs)
            if len(newconfs) > 0:
                self.fill_with_f(self.f, newconfs)
            self._refine_Q()

    def _refine_Q(self):
        oldconfs = copy.copy(self.confs)
        fine = False
        while not fine:
            print("%d th time rotation Q refinement"%(self.Q_refine_count))
            self.Q_refine_count += 1
            fine = True
            for leaf in self.QLeafNodes():
                if leaf.error < self.E_CUTOFF: continue
                node = leaf
                tree = leaf.tree
                # I want to restrict the density in very close
                # if tree.pos[0] < 2.5 and node.size[0] < np.pi/8.0:continue
                node_err = 0
                testgrids = node.testgrid
                #print(node.testgrid)
                for g in testgrids:
                    g_iterp = tree.interpolation(g.q, node=node)
                    err = np.abs(g_iterp - g.values)
                    err = err[0]
                    if err > node_err:
                        node_err = err
                        self.max_err_conf = g
                if node.error == 0:continue
                if node_err < node.error: node.error = node_err
                if node.error > self.E_CUTOFF:
                    printStr=("\nmax  error  is %5.2f\n"%(node.error)+
                               "conf:%15s"%(self.max_err_conf.idx)+
                              " %5.2f"*3%tuple(self.max_err_conf.loc)+
                              " %5.2f"*4%tuple(self.max_err_conf.q) +
                              '\n' +
                              "area of node %6.3f/4*pi\n"%(node.size[0]) +
                              "conf values is " +
                              " %5.2f"*7%tuple(self.max_err_conf.values)
                               ) 
                    #rint(printStr)
                    tree. subdivideNode(node)
                    fine = False 
            self.confs.update(self._iter_conf())
            newconfs = self.confs.difference(oldconfs)
            print("total %15d confs, %15d new confs\n"%(len(self.confs),len(newconfs)))
            oldconfs = copy.copy(self.confs)
            if len(newconfs) > 0:
                self.fill_with_f(self.f, newconfs)
        
        
    def interpolate(self, R, q):
        return self.mesh.interpolation(R, q)

    def load(self, filename):
        load_count = 0
        import pdb
        with open(filename,'r') as f:
            for i in f:
                if load_count%1000 == 0: print('>'*8+"Looded %10d conf"%(load_count)+'<'*8)
                load_count += 1
                if i[0] == '#':continue
                i = i.split()
                if len(i) != 8:continue
                idx = i[0]
                #if idx=="wtr_wtrR360Q60C0":pdb.set_trace()
                values = np.array([float(v) for v in i[1:]])
                self.mesh.fill(idx, values)
    
    def save(self, confs=None, filename='mesh.dat'):
        if confs is None:confs=self._iter_conf()
        if self.database_name != None: filename = self.database_name
        with open(filename, 'a') as f:
            for conf in confs:
                conf_str='%s'%(conf.idx) + ' %f'*7%tuple(conf.values) + '\n'
                f.write(conf_str)

    def gen_x(self):
        for conf in self._iter_conf(self):
            yield (conf.loc, conf.q)
    
    def _iter_conf(self):
        for Qtree in self.mesh.allgrids.values():
            for conf in Qtree.allgrids.values():
                yield conf

    def QLeafNodes(self):
        for RNode in self.RLeafNodes():
            for Qtree in RNode.grids:
                for n in Qtree.leafNodes.values():
                    n.tree = Qtree
                    yield n

    def RLeafNodes(self):
        for n in self.mesh.leafNodes.values():
            n.tree = self.mesh
            yield n
        
            
        
if __name__ == "__main__":
    pass
