import plotly.graph_objects as go
import plotly.io as pio
from copy import deepcopy

template = deepcopy(pio.templates['simple_white'])

template.layout.annotationdefaults.font.size = 18
template.layout.annotationdefaults.font.color = 'black'
# template.layout.annotationdefaults.font.family = 'Times New Roman'
template.layout.annotationdefaults.arrowhead = 1
template.layout.annotationdefaults.arrowwidth = 2
template.layout.annotationdefaults.arrowsize = 0.5

template.layout.legend.font.size = 36
template.layout.legend.bgcolor='rgba(255,255,255,0.8)'
template.layout.legend.borderwidth = 2
template.layout.legend.bordercolor = 'slategray'

template.layout.shapedefaults.line.width = 2
template.layout.shapedefaults.opacity = 1
template.layout.shapedefaults.fillcolor = None

template.layout.font.family = 'Times New Roman'

template.layout.margin.l = 60
template.layout.margin.r = 60
template.layout.margin.t = 60
template.layout.margin.b = 60

template.layout.coloraxis.colorbar.outlinewidth = 0
template.layout.coloraxis.colorbar.ticklen = 0
template.layout.coloraxis.colorbar.tickfont.size = 18
template.layout.coloraxis.colorbar.xpad = 0
template.layout.coloraxis.colorbar.ypad = 0
template.layout.coloraxis.colorbar.title.side = 'right'
template.layout.coloraxis.colorbar.title.font.size = 24

template.layout.xaxis.mirror = True
template.layout.xaxis.showgrid = True
template.layout.xaxis.ticks = 'inside'
template.layout.xaxis.tickwidth = 2
template.layout.xaxis.tickfont.size = 36
template.layout.xaxis.minor.ticks = 'inside'
template.layout.xaxis.minor.ticklen = 3
template.layout.xaxis.minor.tickwidth = 1
template.layout.xaxis.minor.tickcolor = 'slategray'
template.layout.xaxis.linewidth = 2
template.layout.xaxis.title.font.size = 48

template.layout.yaxis.mirror = True
template.layout.yaxis.showgrid = True
template.layout.yaxis.ticks = 'inside'
template.layout.yaxis.tickwidth = 2
template.layout.yaxis.tickfont.size = 36
template.layout.yaxis.minor.ticks = 'inside'
template.layout.yaxis.minor.ticklen = 3
template.layout.yaxis.minor.tickwidth = 1
template.layout.yaxis.minor.tickcolor = 'slategray'
template.layout.yaxis.linewidth = 2
template.layout.yaxis.title.font.size = 48

template.layout.scene.xaxis.mirror = True
template.layout.scene.xaxis.showgrid = True
template.layout.scene.xaxis.ticks = 'inside'
template.layout.scene.xaxis.tickwidth = 4
template.layout.scene.xaxis.ticklen = 10
template.layout.scene.xaxis.tickfont.size = 15
template.layout.scene.xaxis.linewidth = 4
template.layout.scene.xaxis.title.font.size = 24
template.layout.scene.xaxis.backgroundcolor = '#E5ECF6'
template.layout.scene.xaxis.gridcolor = 'black'
template.layout.scene.xaxis.linecolor = 'black'

template.layout.scene.yaxis.mirror = True
template.layout.scene.yaxis.showgrid = True
template.layout.scene.yaxis.ticks = 'inside'
template.layout.scene.yaxis.tickwidth = 4
template.layout.scene.yaxis.ticklen = 10
template.layout.scene.yaxis.tickfont.size = 15
template.layout.scene.yaxis.linewidth = 4
template.layout.scene.yaxis.title.font.size = 24
template.layout.scene.yaxis.backgroundcolor = '#E5ECF6'
template.layout.scene.yaxis.gridcolor = 'black'
template.layout.scene.yaxis.linecolor = 'black'

template.layout.scene.zaxis.mirror = True
template.layout.scene.zaxis.showgrid = True
template.layout.scene.zaxis.ticks = 'inside'
template.layout.scene.zaxis.tickwidth = 4
template.layout.scene.zaxis.ticklen = 10
template.layout.scene.zaxis.tickfont.size = 15
template.layout.scene.zaxis.linewidth = 4
template.layout.scene.zaxis.title.font.size = 24
template.layout.scene.zaxis.backgroundcolor = '#E5ECF6'
template.layout.scene.zaxis.gridcolor = 'black'
template.layout.scene.zaxis.linecolor = 'black'

pio.templates['thesis'] = template
pio.templates.default = 'thesis'