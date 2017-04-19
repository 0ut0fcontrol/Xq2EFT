"""A mesh with refinement proces.
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
##
from basicMesh import mesh
import numpy as np
import pdb
import copy
from time import time
np.seterr(invalid='warn')

class AdaptMesh(mesh):
    def __init__(self):
        mesh.__init__(self)

    def refine(self, f, toPDB, err_cutoff=1.0, filename='mesh.dat'):
        self.confs.update(self._iter_conf())
        self.n = len(self.confs)
        self.database_name = filename # set global filename for self.save
        self.Q_refine_count = 0
        self.R_refine_count = 0
        self.max_err_conf = None
        self.toPDB = toPDB
        self.f = f
        if len(self.confs) < 1000:
            self.fill_with_f(self.f)
        self.E_CUTOFF = err_cutoff
        self._refine_Q()
        self._refine_R()
        self.n = len(self.confs)
        self.save()
        self.database_name = None

    def _refine_R(self):
        oldconfs = copy.copy(self.confs)
        fine = False
        while not fine:
            print("\n%d th time translocation R refinement"%(self.R_refine_count))
            self.R_refine_count += 1
            fine = True
            leaves = set(self.RLeafNodes())
            for leaf in leaves:
                if leaf.error < self.E_CUTOFF: continue
                # I want to restrict the density in very close
                min_size = 0.4
                if leaf.pos[0] < 6.0:min_size = 0.2
                if leaf.pos[0] < 3.5:min_size = 0.0618
                if leaf.size[0] < min_size:
                    print('\n#'*4 + 'R node in %8.5f A Esacape for size < %3.2f, Node.error:%8.3f'%(leaf.pos[0], min_size,leaf.parent.error)+ '#'*4)
                    continue
                    #print("node.idx: %10s pos:%8.5f %8.5f %8.5f"%(leaf.idx,leaf.size[0],leaf.size[1],leaf.size[2]))
                    #continue # 10/(2**6) = 0.156
                node = leaf
                tree = leaf.tree
                node_err = 0
                testgrids = set()
                for Qtree in node.testgrid:
                    for conf in Qtree.allgrids.values():
                        testgrids.add(conf)
                #min_values = 100.0
                #for g in testgrids:
                #    if g.values[0] <  min_values: min_values = g.values[0]
                #if min_values > 25.0 : continue
                #pdb.set_trace()
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
                              "conf:% 15s"%(self.max_err_conf.idx)+
                              " %5.2f"*3%tuple(self.max_err_conf.loc)+
                              " %5.2f"*4%tuple(self.max_err_conf.q) +
                              '\n' +
                              "size of R node %6.3f\n"%(node.size[0]) +
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
            leaves = set(self.QLeafNodes())
            for leaf in leaves:
                if leaf.error < self.E_CUTOFF: continue
                node = leaf
                tree = leaf.tree
                # I want to restrict the density in very close
                #if tree.pos[0] < 2.5 and node.size[0] < np.pi/8.0:continue
                #if leaf.testgrid[0].values[0] > 99: continue
                min_size = np.pi/2.
                if leaf.tree.pos[0] < 6.0:min_size = np.pi/4.
                if leaf.tree.pos[0] < 3.5:min_size = np.pi/16.
                #if leaf.tree.pos[0] < 2.7:min_size = np.pi/16.
                if leaf.testgrid[0].values[0] >  100 and leaf.size[0] < np.pi/2.0: continue
                #if leaf.testgrid[0].values[0] >  25 and leaf.size[0] < np.pi/4.0: continue
                if leaf.size[0] < min_size:
                    print('#'*4 + '\nQ node in %8.5f Esacape for size < np.pi/%d, Node.error:%8.3f'%(
                        leaf.tree.pos[0], int(np.pi/min_size), leaf.parent.error)+ '#'*4)
                    continue
                #if leaf.size[0] < np.pi/16:continue
                #if leaf.testgrid[0].values[0] > 50 and leaf.size[0] < np.pi/2.0: continue
                #if leaf.testgrid[0].values[0] > 25 and  leaf.size[0] < np.pi/4.0: continue
                #if leaf.testgrid[0].values[0] > 25.0:continue
                node_err = 0
                testgrids = node.testgrid
                #pdb.set_trace()
                #print(node.testgrid)
                for g in testgrids:
                    #t0 = time()
                    g_iterp = tree.interpolation(g.q, node=node)
                    #t1 = time()
                    #print("time of interp is %.2fs"%(t1-t0))
                    err = np.abs(g_iterp - g.values)
                    err = err[0]
                    if err > node_err:
                        node_err = err
                        self.max_err_conf = g
                #if node.error == 0:continue
                if node_err < node.error: node.error = node_err
                if node.error > self.E_CUTOFF:
                    printStr=("\nmax  error  is %5.2f\n"%(node.error)+
                                "conf:%15s "%(self.max_err_conf.idx)+
                              " %5.2f"*3%tuple(self.max_err_conf.loc)+
                              " %5.2f"*4%tuple(self.max_err_conf.q) +
                              '\n' +
                              "area of Q node %6.3f/4*pi\n"%(node.size[0]) +
                              "conf values is " +
                              " %5.2f"*7%tuple(self.max_err_conf.values)
                               )   
                    print(printStr)
                    #t0 = time()
                    tree. subdivideNode(node)
                    #t1 = time()
                    #print("time of subdivideNode %s is %.8fs"%(node.idx,t1-t0))
                    fine = False 
            #t0 = time()
            self.confs.update(self._iter_conf())
            #t1 = time()
            #print("time of _iter_conf is %.2fs"%(t1-t0))
            newconfs = self.confs.difference(oldconfs)
            print("total %15d confs, %15d new confs\n"%(len(self.confs),len(newconfs)))
            oldconfs = copy.copy(self.confs)
            if len(newconfs) > 0:
                t0 = time()
                self.fill_with_f(self.f, newconfs)
                t1 = time()
                print("time of fill %d conf is %.2fs"%(len(newconfs), t1-t0))
        
if __name__ == "__main__":
    pass
