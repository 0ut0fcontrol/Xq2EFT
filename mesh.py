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
        self._refine_R()
        self._refine_Q()
        self._refine_S()
        self.n = len(self.confs)
        self.save()
        self.database_name = None

    def _refine_R(self):
        oldconfs = copy.copy(self.confs)
        fine = False
        while not fine:
            #print("\n%d th time translocation R refinement"%(self.R_refine_count))
            self.R_refine_count += 1
            fine = True
            leaves = tuple(self.RLeafNodes())
            parents = [ leaf.parent for leaf in leaves]
            leaves = set(parents)
            for leaf in leaves:
                if leaf.error < self.E_CUTOFF: continue
                # I want to restrict the density in very close
                min_size = 0.05
                if leaf.size < min_size:
                    print('\n#'*4 + 'R node in %8.5f A Esacape for size < %3.2f, Node.error:%8.3f'%(
                            leaf.pos, min_size, leaf.parent.error)+ '#'*4)
                    continue
                tree = leaf.children[0].tree
                node_err = 0.
                print("There are ")
                testgrids = tuple(self._iter_conf(leaf.testgrid))
                print(leaf.testgrid)
                print(" %d conf is Sphere, total %d conf"%(len(testgrids), len(oldconfs)))
                for g in testgrids:
                    g_iterp = tree.interpolation(g.loc, g.q, node=leaf, neighbors=leaf.grids)
                    err = np.abs(g_iterp - g.values)
                    #print('iterp value:%f, values: %f'%(g_iterp[0], g.values[0]))
                    err = err[0]
                    if err > node_err:
                        node_err = err
                        self.max_err_conf = g
                if node_err < leaf.error: leaf.error = node_err
                if leaf.error > self.E_CUTOFF:
                    printStr=("max  error  is %5.2f\n"%(leaf.error)+
                               "c onf:% 15s"%( self.max_err_conf.idx)+
                              " %5.2f"*3%tuple(self.max_err_conf.loc)+
                              " %5.2f"*4%tuple(self.max_err_conf.q) +
                              '\n' +
                              "size of S node %6.3f\n"%(leaf.size) +
                              "conf values is " +
                              " %5.2f"*7%tuple(self.max_err_conf.values)
                               )   

                    print(printStr)
                    for child in leaf.children:
                        tree.subdivideNode(child)
                    fine = False 
            self.confs.update(self._iter_conf())
            newconfs = self.confs.difference(oldconfs)
            print("total %15d confs, %15d new confs\n"%(len(self.confs),len(newconfs)))
            oldconfs = copy.copy(self.confs)
            if len(newconfs) > 0:
                self.fill_with_f(self.f, newconfs)
            #self._refine_Q()

    def _refine_S(self):
        oldconfs = copy.copy(self.confs)
        fine = False
        while not fine:
            print("\n%d th time translocation R refinement"%(self.R_refine_count))
            self.R_refine_count += 1
            fine = True
            leaves = tuple(self.SLeafNodes())
            for leaf in leaves:
                if leaf.error < self.E_CUTOFF: continue
                # I want to restrict the density in very close
                tree = leaf.tree
                tree_grid_max = 100
                if len(tree.gDict) > tree_grid_max:
                    print('\n'+ '#'*4 + 'S node in %8.5f A Esacape for grids num > %d, Node.error:%8.3f'%(
                            tree.r, tree_grid_max, leaf.parent.error)+ '#'*4)
                    continue
                #min_size = np.pi/8. # 3.14/32,
                #if leaf.size[0] < min_size:
                #    print('\n#'*4 + 'S node in %8.5f A Esacape for size < %3.2f, Node.error:%8.3f'%(
                #            tree.r, min_size, leaf.parent.error)+ '#'*4)
                #    continue
                node_err = 0.
                testgrids = self._iter_conf(leaf.testgrid)
                for g in testgrids:
                    g_iterp = tree.interpolation(g.loc, g.q, node=leaf, neighbors=leaf.grids)
                    err = np.abs(g_iterp - g.values)
                    err = err[0]
                    if err > node_err:
                        node_err = err
                        self.max_err_conf = g
                if node_err < leaf.error: leaf.error = node_err
                if leaf.error > self.E_CUTOFF:
                    printStr=("max  error  is %5.2f\n"%(leaf.error)+
                               "c onf:% 15s"%( self.max_err_conf.idx)+
                              " %5.2f"*3%tuple(self.max_err_conf.loc)+
                              " %5.2f"*4%tuple(self.max_err_conf.q) +
                              '\n' +
                              "size of S node %6.3f\n"%(leaf.size[0]) +
                              "conf values is " +
                              " %5.2f"*7%tuple(self.max_err_conf.values)
                               )   
                    print(printStr) 
                    tree.subdivideNode(leaf)
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
            print("\n%d th time translocation R refinement"%(self.R_refine_count))
            self.R_refine_count += 1
            fine = True
            leaves = tuple(self.QLeafNodes())
            for leaf in leaves:
                if leaf.error < self.E_CUTOFF: continue
                # I want to restrict the density in very close
                tree = leaf.tree

                tree_grid_max = 500
                if len(tree.gDict) > tree_grid_max:
                    r = np.linalg.norm(tree.xyz)
                    print('\n'+ '#'*4 + 'Q node in %8.5f A Esacape for grids num > %d, Node.error:%8.3f'%(
                            r, tree_grid_max, leaf.parent.error)+ '#'*4)
                    continue

                min_size = np.pi/8.
                #if leaf.size[0] < min_size:
                #    print('\n#'*4 + 'Q node in %8.5f A Esacape for size < %3.2f, Node.error:%8.3f'%(
                #            r, min_size, leaf.parent.error)+ '#'*4)
                #    continue
                node_err = 0.
                testgrids = self._iter_conf(leaf.testgrid)
                for g in testgrids:
                    g_iterp = tree.interpolation(g.q, node=leaf, neighbors=leaf.grids)
                    err = np.abs(g_iterp - g.values)
                    err = err[0]
                    if err > node_err:
                        node_err = err
                        self.max_err_conf = g
                if node_err < leaf.error: leaf.error = node_err
                if leaf.error > self.E_CUTOFF:
                    printStr=("max  error  is %5.2f\n"%(leaf.error)+
                              "conf:% 15s"%(self.max_err_conf.idx)+
                              " %5.2f"*3%tuple(self.max_err_conf.loc)+
                              " %5.2f"*4%tuple(self.max_err_conf.q) +
                              '\n' +
                              "size of Q node %6.3f\n"%(leaf.size[0]) +
                              "conf values is " +
                              " %5.2f"*7%tuple(self.max_err_conf.values)
                               )    
                    print(printStr) 
                    tree.subdivideNode(leaf)
                    fine = False 
            self.confs.update(self._iter_conf())
            newconfs = self.confs.difference(oldconfs)
            print("total %15d confs, %15d new confs\n"%(len(self.confs),len(newconfs)))
            oldconfs = copy.copy(self.confs)
            if len(newconfs) > 0:
                self.fill_with_f(self.f, newconfs)
        
if __name__ == "__main__":
    pass
