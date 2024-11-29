 ; Declare global variables
globals
[
  ; Membrane parameters
  mthickness
  dlthickness
  ; Others
  stepsize
  gard ;variable to store the previous xcor of a turtle
  mpotential ; membrane potential
  Y; to track ln(Nout/Nin)
  maxY ; to find the max ypcor
  Vm ; Membrane potential
  start_time ; Time at which excitation starts
  ap_start_time ; time at which action potential starts
  set_off ; flag for
  Na_Perm ; tracker to confirm that the excitation button works
  ap_Na_Perm ; tracker to confirm that the action potential procedure works
  tpotential ; threshold potential
  capacitance ; capacitance of the membrane

]

turtles-own         ; Declare turtles-own variables
[
  charge ; -/+ 1
  diffusion_coefficient
  permeability ; General permeability
  stepsize_vd ; additional stepsize due to drift velocity
]


to setup

; Clear all plots and variables
  clear-all

; Chose the values of the constants associated with the membrane
  set mthickness 5    ; Membrane thickness in patches
  set dlthickness 10   ; Thickness of the diffuse layer in patches
  set tpotential -0.0375   ; Threshold potential for an action potential
  set capacitance 2000

; set action potential to be false - this is a flag for the action potential procedure
  set set_off false

; Create the world
  ask patches [
    if (pxcor <= -2) [set pcolor 89.8]                                                 ; The extracellular medium is green
    if (pxcor > (- mthickness / 2) and pxcor < (mthickness / 2)) [set pcolor 28]       ; The membrane is orange
    if (pxcor >= 2) [set pcolor 69.8]                                                  ; The intracellular medium is green
  ]
  draw_scalebar

; Create Chloride ions in the extracellular medium
  create-turtles N_Cl  [
  set shape "circle"
  set color yellow
  set size 1
  set charge -1
  set diffusion_coefficient 2
  set permeability 0.04
  setxy (- abs random-xcor)  random-ycor
  ]

  ; Create Sodium ions in the extracellular medium
  create-turtles N_Na [
  set shape "circle"
  set color blue
  set size 1
  set charge 1
  set diffusion_coefficient 1
  set permeability 0.004
  setxy (- abs random-xcor)  random-ycor
  ]

  ; Create Potassium ions in the extracellular medium
  create-turtles N_K [
  set shape "circle"
  set color green
  set size 1
  set charge 1
  set diffusion_coefficient 2
  set permeability 0.1
  setxy (abs random-xcor)  random-ycor
  ]

; Create Amino Acid ions in the extracellular medium
  create-turtles N_AA [
  set shape "circle"
  set color red
  set size 1
  set charge -10
  set diffusion_coefficient 0.005
  set permeability 0
  setxy (abs random-xcor)  random-ycor
  ]

; Set the initial tick value to 0
  reset-ticks
end

to go
; Normal movement here
  ask turtles [
    set gard xcor ; store the current x-coordinate in the variable "gard"
    move                                        ; Move the ion
    if (gard * xcor < 0) [cross_membrane?]      ; If a move brings an ion across the membrane, decide whether to accept this move
  ]

  ; report permeabilities for the two trackers
  Na_tracker ; checks how many sodium ions have a permeability of 0.4
  ap_Na_tracker ; checks how many sodium ions have a permeability of 1
  action_potential ; sends turtle to action potential procedure to check if the conditions have been met (i.e. if threshold potential has been reached)

  ;sends turtle to revert procedure to check whether or not to reset permeability after excitation
  revert

  tick

; Call Vm calculation
  calculate_Vm
  ; setting mpotential to be Vm
  set mpotential Vm

; Next part is to calculate Y
  ;Counting the potassium ions in the extracellular dT
  let N_Kin count turtles with [color = green and xcor > ((mthickness / 2) + dlthickness)]
  ;Counting the potassium ions in the intracellular dT
  let N_Kout count turtles with [color = green and xcor < ((- mthickness / 2) - dlthickness)]
  ;Calculating ln(N_out/N_out) so we can track the relative concentrations
  set Y ln((N_Kout + 1) / (N_Kin + 1))

end


; Move the ion
to move
 ; Movement due to Brownian motion
    set heading random-float 360
    set stepsize (sqrt (4 * diffusion_coefficient * 1))    ; calculate the step size in number of patch
    fd stepsize
 ; Movement due to drift velocity
  if (xcor < ((mthickness / 2) + dlthickness) and xcor > ((- mthickness / 2) - dlthickness))[ ; if the ion is within the dT region
    set heading 90 ; set its direction to right
    set stepsize_vd (37 * diffusion_coefficient * charge * (- mpotential) / ((2 * dlthickness) + mthickness)) ; calculate the drift velocity of the molecule
    fd stepsize_vd ; move molecule by the distance travelled in one tick
  ]
end

; Decide whether the ion is allowed to cross the membrane
to cross_membrane?

  ; Create new permeability due to pump activity for sodium
  let pass (permeability + Na_pump_activity)

; If a move is not accepted, based on the membrane permeability, return the ion to its previous position
  if (color != blue or (color = blue and xcor > gard)) and random-float 1 > permeability [ ; checks if the ion IS NOT a sodium ion moving from the right to left and uses the standard permability to check
    ;.. if the ion is let through
    set xcor gard ; if it is rejected, return the ion to its previous position
    ]

; Code to include only sodium ions moving from right to left and give them the enhanced permeability due to the pump
  if color = blue and xcor < gard and random-float 1 > pass [ ; checks if the ion IS a sodium ion moving from the right to left and uses the adjusted permability to check
    ;.. if the ion is let through
     set xcor gard ; if it is rejected, return the ion to its previous position
    ]
end

;calculate Vm based on the charge on either side of dT
to calculate_Vm
  ;To find charge in the diffusive layer, we can add the charges up from each ion in the diffusive layer (since they have difference charges)
  ;To find d_T
  let d_T (mthickness + (2 * dlthickness)) ; from the first page of the assignment

  ;Counting each type of ion on the inside of the diffusive layer - this is done individually to account for the difference in charge for amino acids (10)
  let Q_in_K count turtles with [color = green and xcor > 0 and xcor < d_T] ; counts potassium ions in the lefthand dT
  let Q_in_Na count turtles with [color = blue and xcor > 0  and xcor < d_T]; counts sodium ions in the lefthand dT
  let Q_in_Cl count turtles with [color = yellow and xcor > 0  and xcor < d_T]; counts chloride ions in the lefthand dT
  let Q_in_AA count turtles with [color = red and xcor > 0  and xcor < d_T]; counts amino acid ions in the lefthand dT
  let Q_in (Q_in_K + Q_in_Na - Q_in_Cl - (10 * Q_in_AA)) ; sum the ions to find the net charge

  ;Counting each type of ion on the outside of the diffusive layer
  let Q_out_K count turtles with [color = green and xcor < 0 and xcor > (- d_T)]; counts potassium ions in the righthand dT
  let Q_out_Na count turtles with [color = blue and xcor < 0 and xcor > (- d_T)]; counts sodium ions in the righthand dT
  let Q_out_Cl count turtles with [color = yellow and xcor < 0 and xcor > (- d_T)]; counts chloride ions in the righthand dT
  let Q_out_AA count turtles with [color = red and xcor < 0 and xcor > (- d_T)]; counts amino acid ions in the righthand dT
  let Q_out (Q_out_K + Q_out_Na - Q_out_Cl - (10 * Q_out_AA)) ; sum the ions to find the net charge

  set Vm ((Q_in - Q_out)/ (2 * capacitance)) ; uses the calculated charges on either side of dT to find the potential
end

to excite ; excitation procedure
ask turtles [ ; ask each ion the if statement, and if the requirements are met, reset the permeability
    if set_off = false [ ;if the action potential has not been set off (to avoid this button overriding the action potential procedure)
if color = blue and xcor < 0 [ ; for only sodium ions in the extracellular medium
set start_time ticks; the time where the excitation starts
set permeability 0.4; set the permablity to be it's new value, that it will hold for the next 2 ticks
      ]
    ]
  ]
end

to revert
ask turtles [; this is called during the go procedure, which means it will always be checking
if set_off = false [ ; to avoid a reset during the action potential
if ticks > (start_time + 2) and color = blue ; check if the current time is 2 ticks past where excitation began,
 [set permeability 0.004 ; reset value of permeability
      ]
    ]
  ]
end

to action_potential
  if ticks > 3000 and Vm >= tpotential and set_off = false[ ; check if the system has been allowed to fall to resting potential,
    ; then check if the threshold potential has been met and the action potential is not currently occurring
    ask turtles [
      if color = blue [ ; for sodium ions
        set permeability 1 ; set the new permeability
        set ap_start_time ticks; the time where the excitation starts
        set set_off true
      ]
    ]
  ]

   if set_off = true and ticks = (ap_start_time + 100 )[ ; after 100 ticks, change the permeability of potassium and revert the permeability of sodium
      ask turtles [
      if color = green [ ; for potassium ions
        set permeability 1 ; set new permeability
      ]
    ]
    ask turtles [
      if color = blue [ ; for sodium ions
        set permeability 0.004 ] ; revert to old permeability
    ]
  ]

  if set_off = true and ticks = (ap_start_time + 200)[ ; after 200 ticks, revert the permebability of sodium
    set set_off false ; reset the flag
    ask turtles [
      if color = green [ ; for potassium ions
        set permeability 0.1 ] ; revert to old permeability
    ]
  ]
end

; following procedures are mostly admin stuff

to Na_tracker ;tracker to see if the sodium ions are actually being excited
  set Na_Perm count turtles with [color = blue and permeability = 0.4]
end

to ap_Na_tracker ;tracker to see if the sodium ions are working according to the action potential
  set ap_Na_Perm count turtles with [color = blue and permeability = 1]
end

to draw_scalebar
  ask patches with [pxcor >= min-pxcor + 2 and pxcor < (min-pxcor + 2 + 10) and pycor = min-pycor + 2] [set pcolor 104] ; draw a 10 patches-long scale bar in the lower left corner
  ask patch (min-pxcor + 11) (min-pycor + 4) [set plabel-color 104 set plabel "10 patches"]  ; add the label for the scale bar
end

to clear_plot
set-current-plot "Membrane Potential"
clear-plot
end
@#$#@#$#@
GRAPHICS-WINDOW
236
39
749
553
-1
-1
5.0
1
10
1
1
1
0
0
1
1
-50
50
0
100
1
1
1
ticks
30.0

BUTTON
24
22
100
74
SETUP
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
126
22
201
74
RUN
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
224
556
760
702
Concentration profile
x (pxl)
c (ion/pxl)
-50.0
50.0
0.0
10.0
true
false
"" ""
PENS
"N_Cl" 1.0 0 -4079321 true "" "histogram [xcor] of turtles with [color = yellow]"
"N_Na" 1.0 0 -13345367 true "" "histogram [xcor] of turtles with [color = blue]"
"N_K" 1.0 0 -10899396 true "" "histogram [xcor] of turtles with [color = green]"
"N_AA" 1.0 0 -2674135 true "" "histogram [xcor] of turtles with [color = red]"

TEXTBOX
269
10
493
54
Extracellular medium
18
0.0
1

TEXTBOX
534
10
750
54
Intracellular medium
18
0.0
1

SLIDER
26
101
198
134
N_Cl
N_Cl
0
4000
4000.0
10
1
NIL
HORIZONTAL

SLIDER
25
151
197
184
N_Na
N_Na
0
4000
2000.0
10
1
NIL
HORIZONTAL

SLIDER
24
201
196
234
N_K
N_K
0
4000
2000.0
10
1
NIL
HORIZONTAL

SLIDER
23
251
195
284
N_AA
N_AA
0
2000
200.0
10
1
NIL
HORIZONTAL

SLIDER
782
39
1000
72
V_e
V_e
-0.1
0.1
0.0
0.01
1
V
HORIZONTAL

MONITOR
781
90
871
135
N_K monitor
Y
17
1
11

MONITOR
783
231
1003
276
NIL
Vm
17
1
11

PLOT
781
296
1372
595
Membrane Potential
Time (ticks)
Vm (V)
0.0
75.0
-0.1
0.0
true
false
"" ""
PENS
"default" 1.0 0 -5825686 true "" "plot Vm"

SLIDER
23
298
195
331
Na_pump_activity
Na_pump_activity
0
0.996
0.469
0.001
1
NIL
HORIZONTAL

BUTTON
25
354
118
387
Excitation
excite
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
784
185
856
218
CLEAR
clear_plot
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
1029
230
1109
275
NIL
Na_perm
17
1
11

MONITOR
1132
229
1221
274
ap_Na_Perm
ap_Na_Perm
17
1
11

@#$#@#$#@
## WHAT IS IT?

This program simulates the diffusion of chloride ions (represented by yellow circles) in the extracellular medium, which is separated from the intracellular medium by a semi-permeable membrane. 

## HOW TO USE IT

Press the Setup button to create the ions. Press the Go button to let them diffuse throughout the system.

## THINGS TO NOTICE

The permeability of the membrane to this type of ion is initially set to 0.4.

## THINGS TO TRY

Change the value of the permeability. Create new types of ions with different properties.

## EXTENDING THE MODEL

Check the text of the assignment for further ideas on how to extend this model.

## COPYRIGHT NOTICE

Copyright 2023 Cécile Fradin. All rights reserved.

Permission to use, modify or redistribute this model is hereby granted, provided that both of the following requirements are followed:
a) this copyright notice is included.
b) this model will not be redistributed for profit without permission from Cécile Fradin. 
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
