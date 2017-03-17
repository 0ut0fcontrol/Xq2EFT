import cPickle, base64
try:
	from SimpleSession.versions.v65 import beginRestore,\
	    registerAfterModelsCB, reportRestoreError, checkVersion
except ImportError:
	from chimera import UserError
	raise UserError('Cannot open session that was saved in a'
	    ' newer version of Chimera; update your version')
checkVersion([1, 11, 2, 41380])
import chimera
from chimera import replyobj
replyobj.status('Restoring session...', \
    blankAfter=0)
replyobj.status('Beginning session restore...', \
    blankAfter=0, secondary=True)
beginRestore()

def restoreCoreModels():
	from SimpleSession.versions.v65 import init, restoreViewer, \
	     restoreMolecules, restoreColors, restoreSurfaces, \
	     restoreVRML, restorePseudoBondGroups, restoreModelAssociations
	molInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVRFyaWJib25JbnNpZGVDb2xvcnECSwNOfYdVCWJhbGxTY2FsZXEDSwNHP9AAAAAAAAB9h1UJcG9pbnRTaXplcQRLA0c/8AAAAAAAAH2HVQVjb2xvcnEFSwNLAH1xBihLAV1xB0sBYUsCXXEISwJhdYdVCnJpYmJvblR5cGVxCUsDSwB9h1UKc3RpY2tTY2FsZXEKSwNHP/AAAAAAAAB9h1UMbW1DSUZIZWFkZXJzcQtdcQwoTk5OZVUMYXJvbWF0aWNNb2RlcQ1LA0sBfYdVCnZkd0RlbnNpdHlxDksDR0AUAAAAAAAAfYdVBmhpZGRlbnEPSwOJfYdVDWFyb21hdGljQ29sb3JxEEsDTn2HVQ9yaWJib25TbW9vdGhpbmdxEUsDSwB9h1UJYXV0b2NoYWlucRJLA4h9h1UKcGRiVmVyc2lvbnETSwNLAn2HVQhvcHRpb25hbHEUfXEVVQhvcGVuZWRBc3EWiIlLAyhVEG91dGxpZXIuMDAzMS5wZGJxF05OSwF0cRh9cRkoKFUQb3V0bGllci4wNTkyLnBkYnEaTk5LAXRxG11xHEsBYShVEG91dGxpZXIuMTM0Mi5wZGJxHU5OSwF0cR5dcR9LAmF1h4dzVQ9sb3dlckNhc2VDaGFpbnNxIEsDiX2HVQlsaW5lV2lkdGhxIUsDRz/wAAAAAAAAfYdVD3Jlc2lkdWVMYWJlbFBvc3EiSwNLAH2HVQRuYW1lcSNLA1gQAAAAb3V0bGllci4xMzQyLnBkYn1xJChYEAAAAG91dGxpZXIuMDU5Mi5wZGJdcSVLAWFYEAAAAG91dGxpZXIuMDAzMS5wZGJdcSZLAGF1h1UPYXJvbWF0aWNEaXNwbGF5cSdLA4l9h1UPcmliYm9uU3RpZmZuZXNzcShLA0c/6ZmZmZmZmn2HVQpwZGJIZWFkZXJzcSldcSoofXErWAUAAABUSVRMRV1xLFg7AAAAVElUTEUgICAgIFRISVMgUERCIEZJTEUgSVMgR0VORVJBVEVEIEJZICJtb2wybW9sIiAtLWxmemhhby5xLWFzfXEuWAUAAABUSVRMRV1xL1g7AAAAVElUTEUgICAgIFRISVMgUERCIEZJTEUgSVMgR0VORVJBVEVEIEJZICJtb2wybW9sIiAtLWxmemhhby5xMGFzfXExWAUAAABUSVRMRV1xMlg7AAAAVElUTEUgICAgIFRISVMgUERCIEZJTEUgSVMgR0VORVJBVEVEIEJZICJtb2wybW9sIiAtLWxmemhhby5xM2FzZVUDaWRzcTRLA0sCSwCGfXE1KEsBSwCGXXE2SwFhSwBLAIZdcTdLAGF1h1UOc3VyZmFjZU9wYWNpdHlxOEsDR7/wAAAAAAAAfYdVEGFyb21hdGljTGluZVR5cGVxOUsDSwJ9h1UUcmliYm9uSGlkZXNNYWluY2hhaW5xOksDiH2HVQdkaXNwbGF5cTtLA4h9h3Uu'))
	resInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQZpbnNlcnRxAksDVQEgfYdVC2ZpbGxEaXNwbGF5cQNLA4l9h1UEbmFtZXEESwNYAwAAAExJR32HVQVjaGFpbnEFSwNYAQAAAEF9h1UOcmliYm9uRHJhd01vZGVxBksDSwJ9h1UCc3NxB0sDiYmGfYdVCG1vbGVjdWxlcQhLA0sAfXEJKEsBTl1xCksBSwGGcQthhksCTl1xDEsCSwGGcQ1hhnWHVQtyaWJib25Db2xvcnEOSwNLA32HVQVsYWJlbHEPSwNYAAAAAH2HVQpsYWJlbENvbG9ycRBLA0sDfYdVCGZpbGxNb2RlcRFLA0sBfYdVBWlzSGV0cRJLA4l9h1ULbGFiZWxPZmZzZXRxE0sDTn2HVQhwb3NpdGlvbnEUXXEVKEsBSwGGcRZLAUsBhnEXSwFLAYZxGGVVDXJpYmJvbkRpc3BsYXlxGUsDiX2HVQhvcHRpb25hbHEafVUEc3NJZHEbSwNK/////32HdS4='))
	atomInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQdyZXNpZHVlcQJLEksDfXEDKEsETl1xBEsGSwaGcQVhhksFTl1xBksMSwaGcQdhhnWHVQh2ZHdDb2xvcnEISxJOfYdVBG5hbWVxCUsSWAEAAABIfXEKWAEAAABPXXELKEsASwNLBksJSwxLD2Vzh1UDdmR3cQxLEol9h1UOc3VyZmFjZURpc3BsYXlxDUsSiX2HVQVjb2xvcnEOSxJLA31xDyhLBF1xEChLA0sGSw9lSwVdcREoSwRLBUsHSwhLEEsRZXWHVQlpZGF0bVR5cGVxEksSiX2HVQZhbHRMb2NxE0sSVQB9h1UFbGFiZWxxFEsSWAAAAAB9h1UOc3VyZmFjZU9wYWNpdHlxFUsSR7/wAAAAAAAAfYdVB2VsZW1lbnRxFksSSwF9cRdLCF1xGChLAEsDSwZLCUsMSw9lc4dVCmxhYmVsQ29sb3JxGUsSSwN9cRpOXXEbKEsDSwRLBUsGSwdLCEsPSxBLEWVzh1UMc3VyZmFjZUNvbG9ycRxLEksDfXEdTl1xHihLA0sESwVLBksHSwhLD0sQSxFlc4dVD3N1cmZhY2VDYXRlZ29yeXEfSxJYBAAAAG1haW59h1UGcmFkaXVzcSBLEkc/8AAAAAAAAH1xIUc/+AAAAAAAAF1xIihLAEsDSwZLCUsMSw9lc4dVCmNvb3JkSW5kZXhxI11xJChLAEsGhnElSwBLBoZxJksASwaGcSdlVQtsYWJlbE9mZnNldHEoSxJOfYdVEm1pbmltdW1MYWJlbFJhZGl1c3EpSxJHAAAAAAAAAAB9h1UIZHJhd01vZGVxKksSSwJ9h1UIb3B0aW9uYWxxK31xLChVDHNlcmlhbE51bWJlcnEtiIhdcS4oSwFLBoZxL0sBSwaGcTBLAUsGhnExZYdVB2JmYWN0b3JxMoiJSxJHAAAAAAAAAAB9h4dVCW9jY3VwYW5jeXEziIlLEkc/8AAAAAAAAH2Hh3VVB2Rpc3BsYXlxNEsSiH2HdS4='))
	bondInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQVjb2xvcnECSwxLA31xA05dcQQoSwJLA0sESwVLCksLZXOHVQVhdG9tc3EFXXEGKF1xByhLBksHZV1xCChLBksIZV1xCShLCUsKZV1xCihLCUsLZV1xCyhLDEsNZV1xDChLDEsOZV1xDShLD0sQZV1xDihLD0sRZV1xDyhLEksTZV1xEChLEksUZV1xEShLFUsWZV1xEihLFUsXZWVVBWxhYmVscRNLDFgAAAAAfYdVCGhhbGZib25kcRRLDIh9h1UGcmFkaXVzcRVLDEc/yZmZoAAAAH2HVQtsYWJlbE9mZnNldHEWSwxOfYdVCGRyYXdNb2RlcRdLDEsBfYdVCG9wdGlvbmFscRh9VQdkaXNwbGF5cRlLDEsCfYd1Lg=='))
	crdInfo = cPickle.loads(base64.b64decode('gAJ9cQEoSwB9cQIoSwBdcQMoRz/0j1wo9cKPR7/xZFocrAgxR8ADwIMSbpeNh3EERz/mbpeNT987R7/aLQ5WBBiTR8AA87ZFocrBh3EFRz/8NT987ZFoR7/k/fO2RaHLR8AJdLxqfvnbh3EGR0AHj1wo9cKPRz/cWhysCDEnR8ANyLQ5WBBih3EHR0ADJul41P30Rz/rItDlYEGJR8ARlocrAgxKh3EIR0ADl41P3ztkR7/U3S8an753R8ALcKPXCj1xh3EJZVUGYWN0aXZlcQpLAHVLAX1xCyhLAF1xDChHv8Q5WBBiTdNHP7fO2RaHKwJHP+141P3ztkaHcQ1Hv4BiTdLxqfxHP+1P3ztkWh1HP/ZJul41P32HcQ5Hv8xqfvnbItFHv+Jul41P3ztHP/m6XjU/fO6HcQ9HP8OVgQYk3S9HP/ZBiTdLxqhHQAayLQ5WBBmHcRBHv9KfvnbItDlHP/crAgxJul5HP//KwIMSbpiHcRFHP/A5WBBiTdNHP/AxJul41P5HQAUcrAgxJumHcRJlaApLAHVLAn1xEyhLAF1xFChHP/DlYEGJN0xHQBAi0OVgQYlHQBnfO2RaHKyHcRVHP/g9cKPXCj1HQBN3ztkWhytHQBlYEGJN0vKHcRZHP/rQ5WBBiTdHQArlYEGJN0xHQBiyLQ5WBBmHcRdHP/9YEGJN0vJHP/qTdLxqfvpHQBmhysCDEm+HcRhHQAAGJN0vGqBHQAPztkWhysFHQBuJN0vGp/CHcRlHQANgQYk3S8dHP/CDEm6XjVBHQBveNT987ZGHcRplaApLAHV1Lg=='))
	surfInfo = {'category': (0, None, {}), 'probeRadius': (0, None, {}), 'pointSize': (0, None, {}), 'name': [], 'density': (0, None, {}), 'colorMode': (0, None, {}), 'useLighting': (0, None, {}), 'transparencyBlendMode': (0, None, {}), 'molecule': [], 'smoothLines': (0, None, {}), 'lineWidth': (0, None, {}), 'allComponents': (0, None, {}), 'twoSidedLighting': (0, None, {}), 'customVisibility': [], 'drawMode': (0, None, {}), 'display': (0, None, {}), 'customColors': []}
	vrmlInfo = {'subid': (0, None, {}), 'display': (0, None, {}), 'id': (0, None, {}), 'vrmlString': [], 'name': (0, None, {})}
	colors = {'Ru': ((0.141176, 0.560784, 0.560784), 1, u'default'), 'Re': ((0.14902, 0.490196, 0.670588), 1, u'default'), 'Rf': ((0.8, 0, 0.34902), 1, u'default'), 'Ra': ((0, 0.490196, 0), 1, u'default'), 'Rb': ((0.439216, 0.180392, 0.690196), 1, u'default'), 'Rn': ((0.258824, 0.509804, 0.588235), 1, u'default'), 'Rh': ((0.0392157, 0.490196, 0.54902), 1, u'default'), 'Be': ((0.760784, 1, 0), 1, u'default'), 'Ba': ((0, 0.788235, 0), 1, u'default'), 'Bh': ((0.878431, 0, 0.219608), 1, u'default'), 'Bi': ((0.619608, 0.309804, 0.709804), 1, u'default'), 'Bk': ((0.541176, 0.309804, 0.890196), 1, u'default'), 'Br': ((0.65098, 0.160784, 0.160784), 1, u'default'), 'H': ((1, 1, 1), 1, u'default'), 'P': ((1, 0.501961, 0), 1, u'default'), 'Os': ((0.14902, 0.4, 0.588235), 1, u'default'), 'Ge': ((0.4, 0.560784, 0.560784), 1, u'default'), 'Gd': ((0.270588, 1, 0.780392), 1, u'default'), 'Ga': ((0.760784, 0.560784, 0.560784), 1, u'default'), 'Pr': ((0.85098, 1, 0.780392), 1, u'default'), 'Pt': ((0.815686, 0.815686, 0.878431), 1, u'default'), 'Pu': ((0, 0.419608, 1), 1, u'default'),
'C': ((0.564706, 0.564706, 0.564706), 1, u'default'), 'Pb': ((0.341176, 0.34902, 0.380392), 1, u'default'), 'Pa': ((0, 0.631373, 1), 1, u'default'), 'Pd': ((0, 0.411765, 0.521569), 1, u'default'), 'Cd': ((1, 0.85098, 0.560784), 1, u'default'), 'Po': ((0.670588, 0.360784, 0), 1, u'default'), 'Pm': ((0.639216, 1, 0.780392), 1, u'default'), 'Hs': ((0.901961, 0, 0.180392), 1, u'default'), 'Ho': ((0, 1, 0.611765), 1, u'default'), 'Hf': ((0.301961, 0.760784, 1), 1, u'default'), 'Hg': ((0.721569, 0.721569, 0.815686), 1, u'default'), 'He': ((0.85098, 1, 1), 1, u'default'), 'Md': ((0.701961, 0.0509804, 0.65098), 1, u'default'), 'Mg': ((0.541176, 1, 0), 1, u'default'), 'K': ((0.560784, 0.25098, 0.831373), 1, u'default'), 'Mn': ((0.611765, 0.478431, 0.780392), 1, u'default'), 'O': ((1, 0.0509804, 0.0509804), 1, u'default'), 'Mt': ((0.921569, 0, 0.14902), 1, u'default'), 'S': ((1, 1, 0.188235), 1, u'default'), 'W': ((0.129412, 0.580392, 0.839216), 1, u'default'), 'sky blue': ((0.529412, 0.807843, 0.921569), 1, u'default'), 'cornflower blue': ((0.392157, 0.584314, 0.929412), 1, u'default'),
'plum': ((0.866667, 0.627451, 0.866667), 1, u'default'), 'Eu': ((0.380392, 1, 0.780392), 1, u'default'), 'Zr': ((0.580392, 0.878431, 0.878431), 1, u'default'), 'Er': ((0, 0.901961, 0.458824), 1, u'default'), 'Ni': ((0.313725, 0.815686, 0.313725), 1, u'default'), 'No': ((0.741176, 0.0509804, 0.529412), 1, u'default'), 'Na': ((0.670588, 0.360784, 0.94902), 1, u'default'), 'Nb': ((0.45098, 0.760784, 0.788235), 1, u'default'), 'Nd': ((0.780392, 1, 0.780392), 1, u'default'), 'Ne': ((0.701961, 0.890196, 0.960784), 1, u'default'), 'Np': ((0, 0.501961, 1), 1, u'default'), 'Fr': ((0.258824, 0, 0.4), 1, u'default'), 'Fe': ((0.878431, 0.4, 0.2), 1, u'default'), 'Fm': ((0.701961, 0.121569, 0.729412), 1, u'default'), 'B': ((1, 0.709804, 0.709804), 1, u'default'), 'F': ((0.564706, 0.878431, 0.313725), 1, u'default'), 'Sr': ((0, 1, 0), 1, u'default'), 'Zn': ((0.490196, 0.501961, 0.690196), 1, u'default'), 'N': ((0.188235, 0.313725, 0.972549), 1, u'default'), 'Kr': ((0.360784, 0.721569, 0.819608), 1, u'default'), 'Si': ((0.941176, 0.784314, 0.627451), 1, u'default'),
'Sn': ((0.4, 0.501961, 0.501961), 1, u'default'), 'Sm': ((0.560784, 1, 0.780392), 1, u'default'), 'V': ((0.65098, 0.65098, 0.670588), 1, u'default'), 'Sc': ((0.901961, 0.901961, 0.901961), 1, u'default'), 'Sb': ((0.619608, 0.388235, 0.709804), 1, u'default'), 'Sg': ((0.85098, 0, 0.270588), 1, u'default'), 'Se': ((1, 0.631373, 0), 1, u'default'), 'Co': ((0.941176, 0.564706, 0.627451), 1, u'default'), 'Cm': ((0.470588, 0.360784, 0.890196), 1, u'default'), 'Cl': ((0.121569, 0.941176, 0.121569), 1, u'default'), 'Ca': ((0.239216, 1, 0), 1, u'default'), 'Cf': ((0.631373, 0.211765, 0.831373), 1, u'default'), 'Ce': ((1, 1, 0.780392), 1, u'default'), 'Xe': ((0.258824, 0.619608, 0.690196), 1, u'default'), 'Tm': ((0, 0.831373, 0.321569), 1, u'default'), 'Cs': ((0.341176, 0.0901961, 0.560784), 1, u'default'), 'Cr': ((0.541176, 0.6, 0.780392), 1, u'default'), 'Cu': ((0.784314, 0.501961, 0.2), 1, u'default'), 'La': ((0.439216, 0.831373, 1), 1, u'default'), 'Li': ((0.8, 0.501961, 1), 1, u'default'), 'Tl': ((0.65098, 0.329412, 0.301961), 1, u'default'), 'Lu': ((0, 0.670588, 0.141176), 1, u'default'),
'Lr': ((0.780392, 0, 0.4), 1, u'default'), 'Th': ((0, 0.729412, 1), 1, u'default'), 'Ti': ((0.74902, 0.760784, 0.780392), 1, u'default'), 'tan': ((0.823529, 0.705882, 0.54902), 1, u'default'), 'Te': ((0.831373, 0.478431, 0), 1, u'default'), 'Tb': ((0.188235, 1, 0.780392), 1, u'default'), 'Tc': ((0.231373, 0.619608, 0.619608), 1, u'default'), 'Ta': ((0.301961, 0.65098, 1), 1, u'default'), 'Yb': ((0, 0.74902, 0.219608), 1, u'default'), 'Db': ((0.819608, 0, 0.309804), 1, u'default'), 'Dy': ((0.121569, 1, 0.780392), 1, u'default'), 'At': ((0.458824, 0.309804, 0.270588), 1, u'default'), 'I': ((0.580392, 0, 0.580392), 1, u'default'), 'U': ((0, 0.560784, 1), 1, u'default'), 'Y': ((0.580392, 1, 1), 1, u'default'), 'Ac': ((0.439216, 0.670588, 0.980392), 1, u'default'), 'Ag': ((0.752941, 0.752941, 0.752941), 1, u'default'), 'Ir': ((0.0901961, 0.329412, 0.529412), 1, u'default'), 'Am': ((0.329412, 0.360784, 0.94902), 1, u'default'), 'Al': ((0.74902, 0.65098, 0.65098), 1, u'default'), 'As': ((0.741176, 0.501961, 0.890196), 1, u'default'), 'Ar': ((0.501961, 0.819608, 0.890196), 1, u'default'),
'Au': ((1, 0.819608, 0.137255), 1, u'default'), 'Es': ((0.701961, 0.121569, 0.831373), 1, u'default'), 'In': ((0.65098, 0.458824, 0.45098), 1, u'default'), 'Mo': ((0.329412, 0.709804, 0.709804), 1, u'default')}
	materials = {u'default': ((0.85, 0.85, 0.85), 30)}
	pbInfo = {'category': [u'distance monitor'], 'bondInfo': [{'color': (5, None, {}), 'atoms': [[16, 13], [16, 11], [22, 20], [18, 6], [12, 18]], 'label': (5, u'2.335\xc5', {u'0.846\xc5': [0], u'0.797\xc5': [1], u'0.000\xc5': [3], u'1.168\xc5': [2]}), 'halfbond': (5, False, {}), 'labelColor': (5, None, {}), 'labelOffset': (5, chimera.Vector(-1e+99, 0.0, 0.0), {chimera.Vector(-1e+99, 0.0, 0.0): [3], chimera.Vector(-1e+99, 0.0, 0.0): [1], chimera.Vector(-1e+99, 0.0, 0.0): [0], chimera.Vector(-1e+99, 0.0, 0.0): [4]}), 'drawMode': (5, 0, {}), 'display': (5, 2, {})}], 'lineType': (1, 2, {}), 'color': (1, 6, {}), 'optional': {'fixedLabels': (True, False, (1, False, {}))}, 'display': (1, True, {}), 'showStubBonds': (1, False, {}), 'lineWidth': (1, 1, {}), 'stickScale': (1, 1, {}), 'id': [-2]}
	modelAssociations = {}
	colorInfo = (8, (u'H', (1, 1, 1, 1)), {(u'green', (0, 1, 0, 1)): [7], (u'O', (1, 0.0509804, 0.0509804, 1)): [4], (u'sky blue', (0.529412, 0.807843, 0.921569, 1)): [1], (u'tan', (0.823529, 0.705882, 0.54902, 1)): [0], (u'plum', (0.866667, 0.627451, 0.866667, 1)): [2], (u'gray', (0.745, 0.745, 0.745, 1)): [3], (u'yellow', (1, 1, 0, 1)): [6]})
	viewerInfo = {'cameraAttrs': {'center': (-1.3921506511017, 0.021210263279751, 3.5953628871827), 'fieldOfView': 18.537294142553, 'nearFar': (7.8705550310897, -0.071197851695135), 'ortho': False, 'eyeSeparation': 50.8, 'focal': -3.183}, 'viewerAttrs': {'silhouetteColor': None, 'clipping': False, 'showSilhouette': False, 'showShadows': False, 'viewSize': 5.5016586772818, 'labelsOnTop': True, 'depthCueRange': (0.5, 1), 'silhouetteWidth': 2, 'singleLayerTransparency': True, 'shadowTextureSize': 2048, 'backgroundImage': [None, 1, 2, 1, 0, 0], 'backgroundGradient': [('Chimera default', [(1, 1, 1, 1), (0, 0, 1, 1)], 1), 1, 0, 0], 'depthCue': True, 'highlight': 0, 'scaleFactor': 1.5264412841657, 'angleDependentTransparency': True, 'backgroundMethod': 0}, 'viewerHL': 7, 'cameraMode': 'mono', 'detail': 1.5, 'viewerFog': None, 'viewerBG': None}

	replyobj.status("Initializing session restore...", blankAfter=0,
		secondary=True)
	from SimpleSession.versions.v65 import expandSummary
	init(dict(enumerate(expandSummary(colorInfo))))
	replyobj.status("Restoring colors...", blankAfter=0,
		secondary=True)
	restoreColors(colors, materials)
	replyobj.status("Restoring molecules...", blankAfter=0,
		secondary=True)
	restoreMolecules(molInfo, resInfo, atomInfo, bondInfo, crdInfo)
	replyobj.status("Restoring surfaces...", blankAfter=0,
		secondary=True)
	restoreSurfaces(surfInfo)
	replyobj.status("Restoring VRML models...", blankAfter=0,
		secondary=True)
	restoreVRML(vrmlInfo)
	replyobj.status("Restoring pseudobond groups...", blankAfter=0,
		secondary=True)
	restorePseudoBondGroups(pbInfo)
	replyobj.status("Restoring model associations...", blankAfter=0,
		secondary=True)
	restoreModelAssociations(modelAssociations)
	replyobj.status("Restoring camera...", blankAfter=0,
		secondary=True)
	restoreViewer(viewerInfo)

try:
	restoreCoreModels()
except:
	reportRestoreError("Error restoring core models")

	replyobj.status("Restoring extension info...", blankAfter=0,
		secondary=True)


try:
	import StructMeasure
	from StructMeasure.DistMonitor import restoreDistances
	registerAfterModelsCB(restoreDistances, 1)
except:
	reportRestoreError("Error restoring distances in session")


def restoreMidasBase():
	formattedPositions = {}
	import Midas
	Midas.restoreMidasBase(formattedPositions)
try:
	restoreMidasBase()
except:
	reportRestoreError('Error restoring Midas base state')


def restoreMidasText():
	from Midas import midas_text
	midas_text.aliases = {}
	midas_text.userSurfCategories = {}

try:
	restoreMidasText()
except:
	reportRestoreError('Error restoring Midas text state')


def restore_volume_data():
 volume_data_state = \
  {
   'class': 'Volume_Manager_State',
   'data_and_regions_state': [ ],
   'version': 2,
  }
 from VolumeViewer import session
 session.restore_volume_data_state(volume_data_state)

try:
  restore_volume_data()
except:
  reportRestoreError('Error restoring volume data')


def restore_cap_attributes():
 cap_attributes = \
  {
   'cap_attributes': [ ],
   'cap_color': None,
   'cap_offset': 0.01,
   'class': 'Caps_State',
   'default_cap_offset': 0.01,
   'mesh_style': False,
   'shown': True,
   'subdivision_factor': 1.0,
   'version': 1,
  }
 import SurfaceCap.session
 SurfaceCap.session.restore_cap_attributes(cap_attributes)
registerAfterModelsCB(restore_cap_attributes)

geomData = {'AxisManager': {}, 'CentroidManager': {}, 'PlaneManager': {}}

try:
	from StructMeasure.Geometry import geomManager
	geomManager._restoreSession(geomData)
except:
	reportRestoreError("Error restoring geometry objects in session")


def restoreSession_RibbonStyleEditor():
	import SimpleSession
	import RibbonStyleEditor
	userScalings = [('licorice', [[0.35, 0.35], [0.35, 0.35], [0.35, 0.35], [0.35, 0.35, 0.35, 0.35], [0.35, 0.35]])]
	userXSections = []
	userResidueClasses = []
	residueData = [(3, 'Chimera default', 'rounded', u'unknown'), (4, 'Chimera default', 'rounded', u'unknown'), (5, 'Chimera default', 'rounded', u'unknown')]
	flags = RibbonStyleEditor.NucleicDefault1
	SimpleSession.registerAfterModelsCB(RibbonStyleEditor.restoreState,
				(userScalings, userXSections,
				userResidueClasses, residueData, flags))
try:
	restoreSession_RibbonStyleEditor()
except:
	reportRestoreError("Error restoring RibbonStyleEditor state")

trPickle = 'gAJjQW5pbWF0ZS5UcmFuc2l0aW9ucwpUcmFuc2l0aW9ucwpxASmBcQJ9cQMoVQxjdXN0b21fc2NlbmVxBGNBbmltYXRlLlRyYW5zaXRpb24KVHJhbnNpdGlvbgpxBSmBcQZ9cQcoVQZmcmFtZXNxCEsBVQ1kaXNjcmV0ZUZyYW1lcQlLAVUKcHJvcGVydGllc3EKXXELVQNhbGxxDGFVBG5hbWVxDWgEVQRtb2RlcQ5VBmxpbmVhcnEPdWJVCGtleWZyYW1lcRBoBSmBcRF9cRIoaAhLFGgJSwFoCl1xE2gMYWgNaBBoDmgPdWJVBXNjZW5lcRRoBSmBcRV9cRYoaAhLAWgJSwFoCl1xF2gMYWgNaBRoDmgPdWJ1Yi4='
scPickle = 'gAJjQW5pbWF0ZS5TY2VuZXMKU2NlbmVzCnEBKYFxAn1xA1UHbWFwX2lkc3EEfXNiLg=='
kfPickle = 'gAJjQW5pbWF0ZS5LZXlmcmFtZXMKS2V5ZnJhbWVzCnEBKYFxAn1xA1UHZW50cmllc3EEXXEFc2Iu'
def restoreAnimation():
	'A method to unpickle and restore animation objects'
	# Scenes must be unpickled after restoring transitions, because each
	# scene links to a 'scene' transition. Likewise, keyframes must be 
	# unpickled after restoring scenes, because each keyframe links to a scene.
	# The unpickle process is left to the restore* functions, it's 
	# important that it doesn't happen prior to calling those functions.
	import SimpleSession
	from Animate.Session import restoreTransitions
	from Animate.Session import restoreScenes
	from Animate.Session import restoreKeyframes
	SimpleSession.registerAfterModelsCB(restoreTransitions, trPickle)
	SimpleSession.registerAfterModelsCB(restoreScenes, scPickle)
	SimpleSession.registerAfterModelsCB(restoreKeyframes, kfPickle)
try:
	restoreAnimation()
except:
	reportRestoreError('Error in Animate.Session')

def restoreLightController():
	import Lighting
	Lighting._setFromParams({'ratio': 1.25, 'brightness': 1.16, 'material': [30.0, (0.85, 0.85, 0.85), 1.0], 'back': [(0.35740674433659325, 0.6604015517481454, -0.6604015517481455), (1.0, 1.0, 1.0), 0.0], 'mode': 'two-point', 'key': [(-0.35740674433659325, 0.6604015517481454, 0.6604015517481455), (1.0, 1.0, 1.0), 1.0], 'contrast': 0.83, 'fill': [(0.25056280708573153, 0.25056280708573153, 0.9351131265310293), (1.0, 1.0, 1.0), 0.0]})
try:
	restoreLightController()
except:
	reportRestoreError("Error restoring lighting parameters")


def restoreRemainder():
	from SimpleSession.versions.v65 import restoreWindowSize, \
	     restoreOpenStates, restoreSelections, restoreFontInfo, \
	     restoreOpenModelsAttrs, restoreModelClip, restoreSilhouettes

	curSelIds =  []
	savedSels = []
	openModelsAttrs = { 'cofrMethod': 0 }
	from chimera import Point
	openModelsAttrs['cofr'] = Point(-2.32628, -0.573985, 3.59536)
	windowSize = (1864, 872)
	xformMap = {0: (((0.24085045920787, -0.88643853394015, -0.39524395231983), 31.611106929862), (-5.2556544062118, 0.24110915283208, 5.6843595710697), True), 1: (((-0.21724896247622, -0.52207222470775, 0.82476874364376), 155.76893615651), (-1.2842808932621, 1.6062919677941, 4.0957949649428), True), 2: (((0.28162879842346, -0.039178497638975, 0.95872324746056), 66.883398971617), (-0.868129924858, -1.6236021110072, -3.4048091153303), True)}
	fontInfo = {'face': ('Sans Serif', 'Normal', 16)}
	clipPlaneInfo = {}
	silhouettes = {0: True, 1: True, 2: True, 41: True}

	replyobj.status("Restoring window...", blankAfter=0,
		secondary=True)
	restoreWindowSize(windowSize)
	replyobj.status("Restoring open states...", blankAfter=0,
		secondary=True)
	restoreOpenStates(xformMap)
	replyobj.status("Restoring font info...", blankAfter=0,
		secondary=True)
	restoreFontInfo(fontInfo)
	replyobj.status("Restoring selections...", blankAfter=0,
		secondary=True)
	restoreSelections(curSelIds, savedSels)
	replyobj.status("Restoring openModel attributes...", blankAfter=0,
		secondary=True)
	restoreOpenModelsAttrs(openModelsAttrs)
	replyobj.status("Restoring model clipping...", blankAfter=0,
		secondary=True)
	restoreModelClip(clipPlaneInfo)
	replyobj.status("Restoring per-model silhouettes...", blankAfter=0,
		secondary=True)
	restoreSilhouettes(silhouettes)

	replyobj.status("Restoring remaining extension info...", blankAfter=0,
		secondary=True)
try:
	restoreRemainder()
except:
	reportRestoreError("Error restoring post-model state")
from SimpleSession.versions.v65 import makeAfterModelsCBs
makeAfterModelsCBs()

from SimpleSession.versions.v65 import endRestore
replyobj.status('Finishing restore...', blankAfter=0, secondary=True)
endRestore({})
replyobj.status('', secondary=True)
replyobj.status('Restore finished.')

