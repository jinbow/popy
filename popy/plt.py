'''
a collection of commonly used routines

Created on Jan 28, 2013

@author: Jinbo Wang <jinbow@gmail.com>
@organization: Scripps Institute of Oceanography
'''


def image2mask(imgfn,thh=50):
    """load data from image file imgfn,
    use thh to determine the threshold value for masking, 
    below which (darker) will be masked out. """
    from scipy import misc
    ld = misc.imread(imgfn).astype('uint8')
    print ld
    ld[ld<thh] = 1
    ld[ld>=thh] = 0
    ld=ld.astype('bool')
    print ld
    print "read mask data from figure %s, data size ="%imgfn,ld.shape
    print "mask sure masked region is in black color"

    return ld
  
def image4mask(arr, figfn):
    from scipy import misc
    misc.imsave(figfn,arr)
    print "save numpy data to figure "+figfn
    print arr
    return
 
# def image4mask(data, figfn):
#     import Image
#     import numpy as np
#     rescaled = (255.0 / data.max() * (data - data.min())).astype(np.uint8)
#     im=Image.fromarray(rescaled)
#     im.save(figfn)
#     print "save numpy data to figure "+figfn
#     print data
#     return


def savenpeek():
    from pylab import savefig
    import os
    fn = '/tmp/tmp.png'
    savefig(fn,dpi=50)
    os.popen('eog '+fn)
    return

def trimwhite(figfn):
    import os
    print "exec "+'convert -trim '+figfn+' '+figfn
    os.system('convert -trim '+figfn+' '+figfn)
    return

def setlabel(ax,xlab='',ylab='',title=''):
    """
    set xlabel, ylabel, title for ax
    """
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    return

def drawbox(p):
    from pylab import plot,array
    coor = array([[p[0],p[3]],
                  [p[0],p[2]],
                  [p[1],p[2]],
                  [p[1],p[3]],
                  [p[0],p[3]] ])
    plot(coor[:,0],coor[:,1],color='k')
    return

def savefig(i,frst='',figfn='', caption='',scriptname='',
            SectionTitle='', FrameTitle='',ftex='',
            figpath='../figures/', docpath='../documents/',
            savepdf=False, trim=True, comment='',dpi=150,
            append=True,IncludeScript=False,transparent=False,
            wosavefig=False):
    """
    save figure in png format
    write rst file
    """
    import sys,os
    import pylab as plt
    
    if (not append) and os.path.exists(frst):
        os.remove(frst) 
    
    def tofile(filename, lines, header='', method='',match=''):
        "check if a file exists, if not, create and write initial_words to it"
        
        if not os.path.exists(filename):
            f = open(filename,'w')
            f.writelines(header)
            f.close()
            
        "read lines from a file and close it after"
        f = open(filename,'r')
        a = f.readlines()
        f.close()
        
        if method == "insertbefore":
            if match =='':
                print "Error, unable to insert empty line"
                return
            
            for ni,line in zip(range(len(a)),a):
                if match in line:
                    break
                
        f = open(filename,'w')
        f.writelines(a[:ni])
        f.writelines(lines)
        f.writelines(a[ni:])
        f.close()
        
        return

    # setup file names
    sfn = sys.argv[0].split('/')[-1].split('.')[0].replace('_','-')
    print sfn
    
    if frst == '':
        frst = docpath+sfn+'.rst'
        if i==0 and os.path.exists(frst):
            os.remove(frst)
    
    if ftex =='':
        ftex = docpath+sfn+'.beamer.tex'
        
    if figfn =='':
        figpath =figpath+sfn+'/'
                
        try:
            os.mkdir(figpath)
        except OSError:
            pass
        
        figfn = figpath+sfn+'_%03i.png'%i
        
    if SectionTitle == '':
        SectionTitle = sfn

    #  write to rst file
    # if rst file does not exist, create a file and write the following header
    header = ['.. _'+sfn+':\n\n'+SectionTitle.replace('_','-')+ \
                     '\n============================================================\n\n',
                     '.. note:: prepared by '+sys.argv[0]+'\n']
    if IncludeScript:
        header += [ '.. literalinclude:: ../scripts/'+sys.argv[0] ]
 
    lines = ['\n.. figure:: %s\n'%figfn,
             '    :width: 700px\n',
             '    :align: center\n\n',
             '    '+caption+ "\n\n",
             '..\n '+comment+'\n\n',
             ':ref:`Top <'+sfn+'>`\n\n']

    tofile(frst,header=header,lines=lines,method='insertbefore',match="literalinclude")
    
    # write to beamer tex file
    # header
    header = ["\\documentclass{beamer}\n",
              "\\usepackage{graphicx}\n",
              "\\title{"+SectionTitle+"}\n",
              "\\author{Jinbo Wang}\n",
              "\\date{\\today}\n",
              "\\begin{document}\n\n",
              "\\begin{frame}\n\n",
              "\\maketitle\n",
              "\\end{frame}\n\n",
              "\\end{document}"]
    caption = caption.replace(":math:`",'$')
    caption = caption.replace('`','$')
    lines = ['%======================================================\n'
             '\\begin{frame}{'+FrameTitle+'}\n',
             '\\begin{figure}\n',
             '\\centering\n',
             '\\includegraphics[width=4in]{%s}\n'%figfn,
             '\\caption{'+caption+'}\n',
             '\\end{figure}\n',
             '\\end{frame}\n']
    
    tofile(ftex,header=header,lines=lines,method='insertbefore',match="end{document}")
    if not wosavefig:
        plt.savefig(figfn,dpi=dpi,transparent=transparent)
        if trim:
            trimwhite(figfn)
        if savepdf:
            plt.savefig(figfn.replace('png','pdf'),bbox_inches='tight')#,pad_inches=0.5)
    plt.clf()
    return i+1

def yticknum(ax,n):
    import matplotlib.pylab as plt
    ax.yaxis.set_major_locator(plt.MaxNLocator(n))
    return
def xticknum(ax,n):
    import matplotlib.pylab as plt
    ax.xaxis.set_major_locator(plt.MaxNLocator(n))
    return
def noxtick(ax):
    import matplotlib.pylab as plt
    plt.setp(ax.get_xticklabels(),visible=False)
    return
def noytick(ax):
    import matplotlib.pylab as plt
    plt.setp(ax.get_yticklabels(),visible=False)
    return

def columns(num,p,fig):
    """
    num is the number of columns
    p is the margin width, [left,right,bottom,top,margin]
    """
    ax = []
    dx = (1. - p[0] - p[1] -(num-1)*p[-1])/num
    dy = 1. - p[2] - p[3]
    for i in range(num):
        ax.append(fig.add_axes([p[0]+dx*i+p[-1]*i, p[2],dx,  dy]))
        if i>0:
            ax[-1].set_yticks([])

    return ax


def putpathinfo(n,fig):
    """ Insert the file path info in grey font into the bottom left corner in a figure with
    figure handle 'fig'. The input n is the levels of the path. Nothing fancy,
    just remind yourself which file producted the figure.

    For example: putpathinfo(1,fig) will put level1/filename.py into the figure
    with handle 'fig'.
    """
    import os
    import sys
    import pylab
    ax=fig.add_axes([0,0,0.4,0.02])
    ax.set_axis_off()
    txt= reduce(lambda l,r:l+'/'+r,os.getcwd().split('/')[-n:])+sys.argv[0].replace('./','/')
    if (n==0):
        txt =sys.argv[0].replace('./','/')

    ax.text(0,0.07,txt,transform=pylab.gca().transAxes,color='grey')
    return


def cpath(n):
    """Get the path name, goes back n levels
       For example: you are currently in /home/yourname/level1/level2/level3
       cpath(3) will give you /level1/level2/level3/filename
       Strangely, I forgot whether I wrote this, or somebody else did. -- Jinbo Wang
    """
    import os
    import sys
    return reduce(lambda l,r:l+'/'+r, os.getcwd().split('/')[-n:])+'/'+sys.argv[0]

def path_to_poly(path):
    """convert path to polygon"""
    import matplotlib as mpl
    n = len(path)
    return [mpl.path.Path.to_polygons(path[i])[0] for i in range(n) ]

def path_length(path):
    """calculate the length of a path"""
    import numpy as np
    
    poly = path_to_poly(path)
    length =[]
    for p in poly:
        if p.shape[0] == 2 and p.shape[1]!=2:
            p=p.T
        if p.shape[1] !=2:
            print "############   error in polygon data"
            return
    
        dd = np.diff(p,axis=0)
        length.append(np.sqrt(dd[:,0]**2+dd[:,1]**2).sum())
    return np.array(length)

def findcontour(cs):
    """find contour for certain level and return a path
    cs = contour(x,y,z,[level])"""
    #fig = plt.figure(9999)
    #ax = fig.add_subplot(111)
    #cs = ax.contour(x,y,z,[value])
    path = cs.collections[0].get_paths()
    return path

def maskcontour(x,y,z,cs):
    """create a mask for a contour
    outside of the contour is masked
    """
    from matplotlib import path
    import numpy as np
    
    ny,nx = z.shape
    
    pa = findcontour(cs)
    i = path_length(pa)
    
    pa = pa[np.argmax(i)]

    points = np.c_[x.reshape(-1,1),y.reshape(-1,1)]

    mask = pa.contains_points(points).reshape(ny,nx)
    mask = ~mask
    d = np.ma.array(z,mask = mask)
    return d, mask

def maskpoly(x2d,y2d,poly,s=0):
    """
    x2d(ny,nx), y2d(ny,nx), poly(npoints, 2)
    example:
        indian = np.array([ [80,0],[80,217],[113,267],[160,320],[883,320],[883,0],[80,0]])
        x2d,y2d=np.meshgrid(np.arange(2160),np.arange(320))
        ma_ip = maskpoly(x2d,y2d,indian)
    output: masked array (ny,nx) ; True for points inside the polygon
    """
    import numpy as np
    from matplotlib import path
    from mpl_toolkits.basemap import shiftgrid
    ny,nx=x2d.shape
    if s !=0:
        x2ds, tmp = shiftgrid(x2d[1,s], x2d, x2d[1,:], start=False)
        x2d=x2ds
    points=np.c_[x2d.reshape(-1,1),y2d.reshape(-1,1)]
    mask = path.Path(poly).contains_points(points).reshape(ny,nx)
    return mask

def plotpath(path,ax,maxlength=False):
    """convert path to patch and add to ax"""
    import matplotlib as mpl
    import numpy as np
    if maxlength:
        pathlength = path_length(path)
        i = np.argmax(pathlength)
        patch = mpl.patches.PathPatch(path[i],facecolor='orange',lw=2)
        ax.add_patch(patch)
    else:
        n = len(path)
        for i in range(n):
            patch = mpl.patches.PathPatch(path[i],facecolor='orange',lw=2)
            ax.add_patch(patch)
    return

def TS_diagram(ss,ts,ax,dlev=0.1):
    from seawater import csiro as sw
    from numpy import linspace,meshgrid,arange
    from pylab import clabel
    t= linspace(ts.min(), ts.max(), 30)
    s= linspace(ss.min(), ss.max(), 30)
    s2d,t2d = meshgrid(s,t)
    
    #ax.scatter(ss,ts,c=colors, s=size, facecolor=facecolor, edgecolor = 'none', marker = marker)
    h=ax.contour(s2d,t2d,sw.pden(s2d,t2d,s2d*0)-1000,levels=arange(20,30,dlev),colors='k')
    clabel(h,inline=1,fontsize=9,fmt='%3.1f')
    return
