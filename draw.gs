'reinit'
'open model.ctl'
'set background 1'
'c'
'set dbuff on'

t_cur = 1
t_end = 1200

WHILE (t_cur <= t_end)
    'set t ' t_cur
    
    
    'set gxout shaded'
    'd theta_p'
    
    'set gxout contour'
    'set clab off'
    'd qv_p'
    
    'draw title Time: ' t_cur*2 ' sec'
    
    'swap'
    '! sleep 0.01'
    t_cur = t_cur + 1

ENDWHILE
