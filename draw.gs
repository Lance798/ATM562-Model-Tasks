'reinit'
'open model_40m.ctl'
'open model_400m.ctl'

'set background 1'
'c'
'set dbuff on'

v = -8
max_v = 8
step = 0.5
levels = ''
while (v <= max_v)
    levels = levels ' ' v  
    v = v + step
endwhile

t_cur = 1
t_end = 108

WHILE (t_cur <= t_end)
    'set t ' t_cur

    'set vpage 0 11 4.25 8.5'
    'set mproj off'
    'set dfile 1'
    
    
    'set gxout grfill'
    'set clevs' levels
    'd theta_p.1'
    'run cbarn'
    
    'set gxout contour'
    'set clab on'
    'd w.1'
    
    'draw title 40m Res - Time: ' t_cur*100 ' sec'
    'draw xlab Distance (km)'
    'draw ylab Height (km)'

    'set vpage 0 11 0 4.25'
    'set mproj off'
    'set dfile 2'
    
    
    'set gxout grfill'
    'set clevs' levels
    'd theta_p.2'
    'run cbarn'
    
    'set gxout contour'
    'set clab on'
    'd w.2'
    
    'draw title 400m Res - Time: ' t_cur*100 ' sec'
    'draw xlab Distance (km)'
    'draw ylab Height (km)'
    'gxprint img/img_' t_cur*100 '.png'
    'swap'
    
    t_cur = t_cur + 1
ENDWHILE
