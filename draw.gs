'reinit'
'open model.ctl'
'set background 1'
'c'
'set dbuff on'

t_cur = 1
t_end = 108

WHILE (t_cur <= t_end)
    'set t ' t_cur
    
    
    'set gxout grfill'
    v = -8
    max_v = 8
    step = 0.5
    levels = ''

    while (v <= max_v)
    levels = levels ' ' v  
    v = v + step
    endwhile

    'set clevs' levels
    'd theta_p'
    'run cbarn'

    
    'set gxout contour'
    'set clab on'
    'set cint 1'
    'd maskout(qc*1000, qc*1000-0.1)'

    'draw title Time: ' t_cur*100 ' sec'
    
    'swap'
    t_cur = t_cur+1

ENDWHILE
