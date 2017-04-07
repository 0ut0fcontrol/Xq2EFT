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
from basicMesh import mesh
import numpy as np
import pdb
import copy
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
        self.fill_with_f(self.f)
        self.E_CUTOFF = err_cutoff
        self._refine_Q()
        self._refine_R()
        self.n = len(self.confs)
        self.database_name = None

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
        
if __name__ == "__main__":
    pass
