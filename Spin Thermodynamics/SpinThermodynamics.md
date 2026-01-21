This is the supplementary Mathematica notebook for "**Thermodynamics of ideal spin fluids and pseudo-gauge ambiguity**"

## Useful definitions

```wl
In[]:= $Assumptions = {T0 > 0};
```

```wl
In[]:= nullFunc[vT_, v\[Nu]_, v\[Alpha]2_, vw2_, v\[Alpha]w_] = 0;
```

```wl
In[]:= ZeroFunc[T_, A_] := 0
```

![0zlisszfjahdg](img/0zlisszfjahdg.png)

```wl
In[]:= ScaleVars = {v\[Alpha]2 -> \[ScriptE]^2 v\[Alpha]2, vw2 -> \[ScriptE]^2 vw2, v\[Alpha]w -> \[ScriptE]^2 v\[Alpha]w};
```

## Pseudo-gauge transformation results

In this section we collect the results of pseudo-gauge transformations of the currents. 

### Belinfante tensor (even levels)

```wl
In[]:= nontrivialCoeffsEven = {1, 2, 7, 8, 9, 13, 15, 16, 25, 27, 29, 31};
```

##### Term 1

![11r05web5ea6a](img/11r05web5ea6a.png)

##### Term 2

![0234r1kf7glwz](img/0234r1kf7glwz.png)

##### Term 7

![0p31sxbhfjrfc](img/0p31sxbhfjrfc.png)

##### Term 8

![1tk80ygfplzlk](img/1tk80ygfplzlk.png)

##### Term 9

![00s0fy59ktclk](img/00s0fy59ktclk.png)

##### Term 13

![1g64l3yv3gom7](img/1g64l3yv3gom7.png)

##### Term 15

![0w1wxj8iat0sx](img/0w1wxj8iat0sx.png)

##### Term 16

![0msgkgw4557vf](img/0msgkgw4557vf.png)

##### Term 25

![08fjglspvnd15](img/08fjglspvnd15.png)

##### Term 27

![0w48noj512qaq](img/0w48noj512qaq.png)

##### Term 29

![0wo7ehe255mtc](img/0wo7ehe255mtc.png)

##### Term 31

![170zhv58e0d9b](img/170zhv58e0d9b.png)

##### Final

![06olv6pf8efz1](img/06olv6pf8efz1.png)

### Belinfante tensor (odd levels)

```wl
In[]:= nontrivialCoeffsOdd = {4, 17, 19, 21, 23, 33, 34, 36};
```

##### Term 4

![05rgz4r2ayh01](img/05rgz4r2ayh01.png)

##### Term 17

![0pyeotzxwflxp](img/0pyeotzxwflxp.png)

##### Term 19

![1n5mwhr4w6jpm](img/1n5mwhr4w6jpm.png)

##### Term 21

![1a5caciwpt8l4](img/1a5caciwpt8l4.png)

```wl
In[]:= \[Omega]\[Mu] + O[\[ScriptE]]^2
```

![060vplnj17t41](img/060vplnj17t41.png)

##### Term 23

![00yrlxwnvg9bm](img/00yrlxwnvg9bm.png)

##### Term 33

![1umxmgjdljf7n](img/1umxmgjdljf7n.png)

##### Term 34

![00abnjaipz5mo](img/00abnjaipz5mo.png)

##### Term 36

![1o08h6slpl144](img/1o08h6slpl144.png)

##### Final

![1agd6s6q1dyt1](img/1agd6s6q1dyt1.png)

### Charge current (odd levels)

![0axy5fndu5la7](img/0axy5fndu5la7.png)

![0f6lko0lxqz6x](img/0f6lko0lxqz6x.png)

### Charge current (even levels)

![1mlzge22breiq](img/1mlzge22breiq.png)

![1i4avyg2rwk5w](img/1i4avyg2rwk5w.png)

## Pseudo-gauge transformation derivation (DO NOT COMPILE unless needed)

In this section we present the derivation/verification of pseudo-gauge transformations of the currents. 
You can skip this section and move directly to later ones for efficiency.

### Preparation

```wl
In[]:= ExpansionStrength = 9;
```

![11wwgosmkli8n](img/11wwgosmkli8n.png)

![0gn4796mhuy7v](img/0gn4796mhuy7v.png)

![1erel8ksvzttk](img/1erel8ksvzttk.png)

```wl
In[]:= a2 = a\[Mu] . \[Eta] . a\[Mu];
 \[Omega]2 = \[Omega]\[Mu] . \[Eta] . \[Omega]\[Mu];
 a\[Omega] = a\[Mu] . \[Eta] . \[Omega]\[Mu];
```

```wl
In[]:= vars = {T0, \[Nu], \[Alpha]2, w2, \[Alpha]w};
 VARS = {T, \[Nu], a2/T^2, \[Omega]2/T^2, a\[Omega]/T^2};
```

```wl
In[]:= T\[Mu]\[Nu]ConstructorEven[\[ScriptCapitalE]_, \[ScriptCapitalP]_, \[ScriptCapitalX]aa_, \[ScriptCapitalX]\[Omega]\[Omega]_, \[ScriptCapitalX]a\[Omega]_, \[ScriptCapitalQ]\[ScriptL]_] := (\[ScriptCapitalE] @@ VARS) u\[Mu]\[TensorProduct]u\[Mu] + (\[ScriptCapitalP] @@ VARS) \[CapitalDelta]\[Mu]\[Nu] + (\[ScriptCapitalX]aa @@ VARS) a\[Mu]\[TensorProduct]a\[Mu] + (\[ScriptCapitalX]\[Omega]\[Omega] @@ VARS) \[Omega]\[Mu]\[TensorProduct]\[Omega]\[Mu] + (\[ScriptCapitalX]a\[Omega] @@ VARS ) (\[Omega]\[Mu]\[TensorProduct]a\[Mu] + a\[Mu]\[TensorProduct]\[Omega]\[Mu]) + (\[ScriptCapitalQ]\[ScriptL] @@ VARS) (\[ScriptL]\[Mu]\[TensorProduct]u\[Mu] + u\[Mu]\[TensorProduct]\[ScriptL]\[Mu]);
 T\[Mu]\[Nu]ConstructorOdd[\[ScriptCapitalQ]a_, \[ScriptCapitalQ]\[Omega]_, \[ScriptCapitalX]a\[ScriptL]_, \[ScriptCapitalX]\[Omega]\[ScriptL]_] := (\[ScriptCapitalQ]a @@ VARS) (u\[Mu]\[TensorProduct]a\[Mu] + a\[Mu]\[TensorProduct]u\[Mu]) + (\[ScriptCapitalQ]\[Omega] @@ VARS) (u\[Mu]\[TensorProduct]\[Omega]\[Mu] + \[Omega]\[Mu]\[TensorProduct]u\[Mu]) + (\[ScriptCapitalX]a\[ScriptL] @@ VARS) (\[ScriptL]\[Mu]\[TensorProduct]a\[Mu] + a\[Mu]\[TensorProduct]\[ScriptL]\[Mu]) + (\[ScriptCapitalX]\[Omega]\[ScriptL] @@ VARS) (\[ScriptL]\[Mu]\[TensorProduct]\[Omega]\[Mu] + \[Omega]\[Mu]\[TensorProduct]\[ScriptL]\[Mu])
 J\[Mu]Constructor[\[ScriptCapitalN]_, \[ScriptCapitalJ]a_, \[ScriptCapitalJ]\[Omega]_, \[ScriptCapitalJ]\[ScriptL]_] := (\[ScriptCapitalN] @@ VARS) u\[Mu] + (\[ScriptCapitalJ]a @@ VARS) a\[Mu] + (\[ScriptCapitalJ]\[Omega] @@ VARS) \[Omega]\[Mu] + (\[ScriptCapitalJ]\[ScriptL] @@ VARS) \[ScriptL]\[Mu];
```

![10nl3ye2c2673](img/10nl3ye2c2673.png)

![14jvkktryzmwz](img/14jvkktryzmwz.png)

![1pun0vvmzc1q5](img/1pun0vvmzc1q5.png)

### Charge current

```wl
In[]:= fJ[1] = u\[Mu]\[TensorProduct]a\[Mu] - a\[Mu]\[TensorProduct]u\[Mu];
 fJ[2] = bar\[Mu];
 fJ[3] = u\[Mu]\[TensorProduct]\[ScriptL]\[Mu] - \[ScriptL]\[Mu]\[TensorProduct]u\[Mu];
 fJ[4] = a\[Omega] (u\[Mu]\[TensorProduct]\[Omega]\[Mu] - \[Omega]\[Mu]\[TensorProduct]u\[Mu]);
 fJ[5] = a\[Mu]\[TensorProduct]\[ScriptL]\[Mu] - \[ScriptL]\[Mu]\[TensorProduct]a\[Mu];
 fJ[6] = a\[Omega] (a\[Mu]\[TensorProduct]\[Omega]\[Mu] - \[Omega]\[Mu]\[TensorProduct]a\[Mu]);
```

![00cy8r3mamegn](img/00cy8r3mamegn.png)

```wl
In[]:= J\[Mu]CorrectionAnsatz = M\[Mu]\[Nu] // partial[#] & // Tr[#, Plus, 2] &;
```

```wl
In[]:= \[Delta]J\[Mu]X = J\[Mu]Constructor[\[Delta]\[ScriptCapitalN]X, \[Delta]\[ScriptCapitalJ]aX, \[Delta]\[ScriptCapitalJ]\[Omega]X, \[Delta]\[ScriptCapitalJ]\[ScriptL]X];
 J\[Mu]CorrectionAnsatz - \[Delta]J\[Mu]X + O[\[ScriptE]]^9 // Simplify
```

![0i4eu5vmn1he9](img/0i4eu5vmn1he9.png)

### Belinfante tensor -- independent structures

![1t7rfwj8ljpzw](img/1t7rfwj8ljpzw.png)

```wl
In[]:= fS[3] = u\[Mu]\[TensorProduct]bar\[Mu]\[TensorProduct]u\[Mu] // Symmetriser;
 fS[4] = (u\[Mu]\[TensorProduct]\[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]a\[Mu] + a\[Mu]\[TensorProduct]\[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]u\[Mu]) // Symmetriser;
 fS[5] = (u\[Mu]\[TensorProduct]\[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]a\[Mu] - a\[Mu]\[TensorProduct]\[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]u\[Mu]) // Symmetriser;
 fS[6] = \[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]bar\[Mu] // Transpose[#, {1, 3, 2, 4}] & // Symmetriser;
```

![0hzutzgumey3l](img/0hzutzgumey3l.png)

```wl
In[]:= fS[17] = (u\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]\[ScriptL]\[Mu]\[TensorProduct]u\[Mu] + u\[Mu]\[TensorProduct]\[ScriptL]\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]u\[Mu]) // Symmetriser;
 fS[18] = (u\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]\[ScriptL]\[Mu]\[TensorProduct]u\[Mu] - u\[Mu]\[TensorProduct]\[ScriptL]\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]u\[Mu]) // Symmetriser;
 fS[19] = a\[Omega] (u\[Mu]\[TensorProduct]\[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]\[Omega]\[Mu] + \[Omega]\[Mu]\[TensorProduct]\[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]u\[Mu]) // Symmetriser;
 fS[20] = a\[Omega] (u\[Mu]\[TensorProduct]\[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]\[Omega]\[Mu] - \[Omega]\[Mu]\[TensorProduct]\[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]u\[Mu]) // Symmetriser;
 fS[21] = (u\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]a\[Mu] + a\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]u\[Mu]) // Symmetriser;
 fS[22] = (u\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]a\[Mu] - a\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]u\[Mu]) // Symmetriser;
 fS[23] = (\[ScriptL]\[Mu]\[TensorProduct]\[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]a\[Mu] + a\[Mu]\[TensorProduct]\[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]\[ScriptL]\[Mu]) // Symmetriser;
 fS[24] = (\[ScriptL]\[Mu]\[TensorProduct]\[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]a\[Mu] - a\[Mu]\[TensorProduct]\[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]\[ScriptL]\[Mu]) // Symmetriser;
```

```wl
In[]:= fS[25] = a\[Omega] (u\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]u\[Mu] + u\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]u\[Mu]) // Symmetriser;
 fS[26] = a\[Omega] (u\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]u\[Mu] - u\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]u\[Mu]) // Symmetriser;
 fS[27] = (u\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]\[ScriptL]\[Mu] + \[ScriptL]\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]u\[Mu]) // Symmetriser;
 fS[28] = (u\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]\[ScriptL]\[Mu] - \[ScriptL]\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]u\[Mu]) // Symmetriser;
 fS[29] = (u\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]\[ScriptL]\[Mu] + \[ScriptL]\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]u\[Mu]) // Symmetriser;
 fS[30] = (u\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]\[ScriptL]\[Mu] - \[ScriptL]\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]u\[Mu]) // Symmetriser;
 fS[31] = a\[Omega] (a\[Mu]\[TensorProduct]\[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]\[Omega]\[Mu] + \[Omega]\[Mu]\[TensorProduct]\[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]a\[Mu]) // Symmetriser;
 fS[32] = a\[Omega] (a\[Mu]\[TensorProduct]\[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]\[Omega]\[Mu] - \[Omega]\[Mu]\[TensorProduct]\[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]a\[Mu]) // Symmetriser;
```

```wl
In[]:= fS[33] = a\[Omega] (u\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]\[ScriptL]\[Mu]\[TensorProduct]u\[Mu] + u\[Mu]\[TensorProduct]\[ScriptL]\[Mu]\[TensorProduct]\[Omega]\[Mu]\[TensorProduct]u\[Mu]) // Symmetriser;
 fS[34] = a\[Omega] (u\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]\[Omega]\[Mu] + \[Omega]\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]u\[Mu]) // Symmetriser;
 fS[35] = a\[Omega] (u\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]\[Omega]\[Mu] - \[Omega]\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]a\[Mu]\[TensorProduct]u\[Mu]) // Symmetriser;
 fS[36] = a\[Omega] (\[ScriptL]\[Mu]\[TensorProduct]\[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]\[Omega]\[Mu] + \[Omega]\[Mu]\[TensorProduct]\[CapitalDelta]\[Mu]\[Nu]\[TensorProduct]\[ScriptL]\[Mu]) // Symmetriser;
```

### Belinfante tensor (even levels)

##### Term 1

![1bq2op4bmwklf](img/1bq2op4bmwklf.png)

```wl
In[]:= XX = ImprovementT\[Mu]\[Nu][1] - T\[Mu]\[Nu]ConstructorEven[\[ScriptCapitalE]X[1], \[ScriptCapitalP]X[1], \[ScriptCapitalX]aaX[1], \[ScriptCapitalX]\[Omega]\[Omega]X[1], \[ScriptCapitalX]a\[Omega]X[1], \[ScriptCapitalQ]\[ScriptL]X[1]] + O[\[ScriptE]]^7 /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify;
 % // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

##### Term 2

![1jz82ly3c3uv0](img/1jz82ly3c3uv0.png)

```wl
In[]:= XX = ImprovementT\[Mu]\[Nu][2] - T\[Mu]\[Nu]ConstructorEven[\[ScriptCapitalE]X[2], \[ScriptCapitalP]X[2], \[ScriptCapitalX]aaX[2], \[ScriptCapitalX]\[Omega]\[Omega]X[2], \[ScriptCapitalX]a\[Omega]X[2], \[ScriptCapitalQ]\[ScriptL]X[2]] + O[\[ScriptE]]^7 /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify;
 % // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

##### Term 7

![1rznpy7xbjt2s](img/1rznpy7xbjt2s.png)

```wl
In[]:= ImprovementT\[Mu]\[Nu][7] - T\[Mu]\[Nu]ConstructorEven[\[ScriptCapitalE]X[7], \[ScriptCapitalP]X[7], \[ScriptCapitalX]aaX[7], \[ScriptCapitalX]\[Omega]\[Omega]X[7], \[ScriptCapitalX]a\[Omega]X[7], \[ScriptCapitalQ]\[ScriptL]X[7]] + O[\[ScriptE]]^9 /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

##### Term 8

![1au9eleni5mck](img/1au9eleni5mck.png)

```wl
In[]:= ImprovementT\[Mu]\[Nu][8] - T\[Mu]\[Nu]ConstructorEven[\[ScriptCapitalE]X[8], \[ScriptCapitalP]X[8], \[ScriptCapitalX]aaX[8], \[ScriptCapitalX]\[Omega]\[Omega]X[8], \[ScriptCapitalX]a\[Omega]X[8], \[ScriptCapitalQ]\[ScriptL]X[8]] + O[\[ScriptE]]^9 /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

##### Term 9

![05nj1s5krr7ok](img/05nj1s5krr7ok.png)

```wl
In[]:= ImprovementT\[Mu]\[Nu][9] - T\[Mu]\[Nu]ConstructorEven[\[ScriptCapitalE]X[9], \[ScriptCapitalP]X[9], \[ScriptCapitalX]aaX[9], \[ScriptCapitalX]\[Omega]\[Omega]X[9], \[ScriptCapitalX]a\[Omega]X[9], \[ScriptCapitalQ]\[ScriptL]X[9]] + O[\[ScriptE]]^9 /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

##### Term 13

![0l5fix749x18l](img/0l5fix749x18l.png)

```wl
In[]:= ImprovementT\[Mu]\[Nu][13] - T\[Mu]\[Nu]ConstructorEven[\[ScriptCapitalE]X[13], \[ScriptCapitalP]X[13], \[ScriptCapitalX]aaX[13], \[ScriptCapitalX]\[Omega]\[Omega]X[13], \[ScriptCapitalX]a\[Omega]X[13], \[ScriptCapitalQ]\[ScriptL]X[13]] + O[\[ScriptE]]^9 /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

##### Term 15

![02wwgmk8jbmpv](img/02wwgmk8jbmpv.png)

```wl
In[]:= ImprovementT\[Mu]\[Nu][15] - T\[Mu]\[Nu]ConstructorEven[\[ScriptCapitalE]X[15], \[ScriptCapitalP]X[15], \[ScriptCapitalX]aaX[15], \[ScriptCapitalX]\[Omega]\[Omega]X[15], \[ScriptCapitalX]a\[Omega]X[15], \[ScriptCapitalQ]\[ScriptL]X[15]] /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

##### Term 16

![0yubjc3eipmy2](img/0yubjc3eipmy2.png)

```wl
In[]:= ImprovementT\[Mu]\[Nu][16] - T\[Mu]\[Nu]ConstructorEven[\[ScriptCapitalE]X[16], \[ScriptCapitalP]X[16], \[ScriptCapitalX]aaX[16], \[ScriptCapitalX]\[Omega]\[Omega]X[16], \[ScriptCapitalX]a\[Omega]X[16], \[ScriptCapitalQ]\[ScriptL]X[16]] /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

##### Term 25 (Scaled)

![0nubuvg3mfneh](img/0nubuvg3mfneh.png)

```wl
In[]:= ImprovementT\[Mu]\[Nu]Scaled[25] - T\[Mu]\[Nu]ConstructorEven[\[ScriptCapitalE]X[25] // scaleFunc, \[ScriptCapitalP]X[25] // scaleFunc, \[ScriptCapitalX]aaX[25] // scaleFunc, \[ScriptCapitalX]\[Omega]\[Omega]X[25] // scaleFunc, \[ScriptCapitalX]a\[Omega]X[25] // scaleFunc, \[ScriptCapitalQ]\[ScriptL]X[25] // scaleFunc] + O[\[ScriptE]]^9 /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

##### Term 27

![0qe2ts547uaic](img/0qe2ts547uaic.png)

```wl
In[]:= XX = ImprovementT\[Mu]\[Nu][27] - T\[Mu]\[Nu]ConstructorEven[\[ScriptCapitalE]X[27], \[ScriptCapitalP]X[27], \[ScriptCapitalX]aaX[27], \[ScriptCapitalX]\[Omega]\[Omega]X[27], \[ScriptCapitalX]a\[Omega]X[27], \[ScriptCapitalQ]\[ScriptL]X[27]] + O[\[ScriptE]]^9 /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify;
 % // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

##### Term 29

![1cnjv4g6mv3hd](img/1cnjv4g6mv3hd.png)

```wl
In[]:= ImprovementT\[Mu]\[Nu][29] - T\[Mu]\[Nu]ConstructorEven[\[ScriptCapitalE]X[29], \[ScriptCapitalP]X[29], \[ScriptCapitalX]aaX[29], \[ScriptCapitalX]\[Omega]\[Omega]X[29], \[ScriptCapitalX]a\[Omega]X[29], \[ScriptCapitalQ]\[ScriptL]X[29]] + O[\[ScriptE]]^9 /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify;
 % // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

##### Term 31 (Scaled)

![0bp5tcp104p91](img/0bp5tcp104p91.png)

```wl
In[]:= ImprovementT\[Mu]\[Nu]Scaled[31] - T\[Mu]\[Nu]ConstructorEven[\[ScriptCapitalE]X[31] // scaleFunc, \[ScriptCapitalP]X[31] // scaleFunc, \[ScriptCapitalX]aaX[31] // scaleFunc, \[ScriptCapitalX]\[Omega]\[Omega]X[31] // scaleFunc, \[ScriptCapitalX]a\[Omega]X[31] // scaleFunc, \[ScriptCapitalQ]\[ScriptL]X[31] // scaleFunc] + O[\[ScriptE]]^9 /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

### Belinfante tensor (odd levels)

##### Term 4

![0yequ69gs1x6e](img/0yequ69gs1x6e.png)

```wl
In[]:= ImprovementT\[Mu]\[Nu][4] - T\[Mu]\[Nu]ConstructorOdd[\[ScriptCapitalQ]aX[4], \[ScriptCapitalQ]\[Omega]X[4], \[ScriptCapitalX]a\[ScriptL]X[4], \[ScriptCapitalX]\[Omega]\[ScriptL]X[4]] + O[\[ScriptE]]^8 /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

##### Term 17

![083zjfk551qkl](img/083zjfk551qkl.png)

```wl
In[]:= ImprovementT\[Mu]\[Nu][17] - T\[Mu]\[Nu]ConstructorOdd[\[ScriptCapitalQ]aX[17], \[ScriptCapitalQ]\[Omega]X[17], \[ScriptCapitalX]a\[ScriptL]X[17], \[ScriptCapitalX]\[Omega]\[ScriptL]X[17]] + O[\[ScriptE]]^10 /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

##### Term 19

![10bbv6vzqtou9](img/10bbv6vzqtou9.png)

```wl
In[]:= ImprovementT\[Mu]\[Nu]Scaled[19] - T\[Mu]\[Nu]ConstructorOdd[\[ScriptCapitalQ]aX[19] // scaleFunc, \[ScriptCapitalQ]\[Omega]X[19] // scaleFunc, \[ScriptCapitalX]a\[ScriptL]X[19] // scaleFunc, \[ScriptCapitalX]\[Omega]\[ScriptL]X[19] // scaleFunc] + O[\[ScriptE]]^8 /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify;
 % // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

##### Term 21

![06jkp91dptkpf](img/06jkp91dptkpf.png)

```wl
In[]:= ImprovementT\[Mu]\[Nu][21] - T\[Mu]\[Nu]ConstructorOdd[\[ScriptCapitalQ]aX[21], \[ScriptCapitalQ]\[Omega]X[21], \[ScriptCapitalX]a\[ScriptL]X[21], \[ScriptCapitalX]\[Omega]\[ScriptL]X[21]] + O[\[ScriptE]]^10 /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify;
 % // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

##### Term 23

![0w198xxrg1xjh](img/0w198xxrg1xjh.png)

```wl
In[]:= ImprovementT\[Mu]\[Nu][23] - T\[Mu]\[Nu]ConstructorOdd[\[ScriptCapitalQ]aX[23], \[ScriptCapitalQ]\[Omega]X[23], \[ScriptCapitalX]a\[ScriptL]X[23], \[ScriptCapitalX]\[Omega]\[ScriptL]X[23]] + O[\[ScriptE]]^10 /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify;
 % // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

##### Term 33

![1d4ni5588ohef](img/1d4ni5588ohef.png)

```wl
In[]:= ImprovementT\[Mu]\[Nu]Scaled[33] - T\[Mu]\[Nu]ConstructorOdd[\[ScriptCapitalQ]aX[33] // scaleFunc, \[ScriptCapitalQ]\[Omega]X[33] // scaleFunc, \[ScriptCapitalX]a\[ScriptL]X[33] // scaleFunc, \[ScriptCapitalX]\[Omega]\[ScriptL]X[33] // scaleFunc] + O[\[ScriptE]]^8 /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify;
 % // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

##### Term 34

![1trkd4aaktpzi](img/1trkd4aaktpzi.png)

```wl
In[]:= ImprovementT\[Mu]\[Nu]Scaled[34] - T\[Mu]\[Nu]ConstructorOdd[\[ScriptCapitalQ]aX[34] // scaleFunc, \[ScriptCapitalQ]\[Omega]X[34] // scaleFunc, \[ScriptCapitalX]a\[ScriptL]X[34] // scaleFunc, \[ScriptCapitalX]\[Omega]\[ScriptL]X[34] // scaleFunc] + O[\[ScriptE]]^8 /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify;
 % // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

##### Term 36

![1u2xisccd52dn](img/1u2xisccd52dn.png)

```wl
In[]:= ImprovementT\[Mu]\[Nu]Scaled[36] - T\[Mu]\[Nu]ConstructorOdd[\[ScriptCapitalQ]aX[36] // scaleFunc, \[ScriptCapitalQ]\[Omega]X[36] // scaleFunc, \[ScriptCapitalX]a\[ScriptL]X[36] // scaleFunc, \[ScriptCapitalX]\[Omega]\[ScriptL]X[36] // scaleFunc] + O[\[ScriptE]]^8 /. \[CapitalLambda] -> \[CapitalLambda]F // Simplify;
 % // MatrixForm
```

|  |  |  |  |
| - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

## Thermodynamics

In this section we derive various core results of the paper concerning thermodynamic constraints and thermodynamic relations in conserved currents.

### Thermodynamic Constraints

![1bxb75z8q8fdt](img/1bxb75z8q8fdt.png)

![07bfw5mqugbg6](img/07bfw5mqugbg6.png)

![143mzdtsr28tt](img/143mzdtsr28tt.png)

```wl
In[]:= \[Delta]Consti = {\[Delta]\[ScriptCapitalE], \[Delta]\[ScriptCapitalP], \[Delta]\[ScriptCapitalX]aa, \[Delta]\[ScriptCapitalX]\[Omega]\[Omega], \[Delta]\[ScriptCapitalX]a\[Omega], \[Delta]\[ScriptCapitalQ]\[ScriptL], \[Delta]\[ScriptCapitalN], \[Delta]\[ScriptCapitalN]\[ScriptL]};
```

### Conservation Identities

![0w6jh5qk4yvaw](img/0w6jh5qk4yvaw.png)

```wl
Out[]= {0, 0}
```

```wl
In[]:= TraceIdentity[\[ScriptCapitalE]_, \[ScriptCapitalP]_, \[ScriptCapitalX]aa_, \[ScriptCapitalX]\[Omega]\[Omega]_, \[ScriptCapitalX]a\[Omega]_, \[ScriptCapitalQ]\[ScriptL]_, \[ScriptCapitalN]_, \[ScriptCapitalN]\[ScriptL]_] := -\[ScriptCapitalE] + 3 \[ScriptCapitalP] + \[ScriptCapitalX]aa T0^2 v\[Alpha]2 + \[ScriptCapitalX]\[Omega]\[Omega] T0^2 vw2 + 2 \[ScriptCapitalX]a\[Omega] T0^2 v\[Alpha]w
```

![19t66h64e9wrc](img/19t66h64e9wrc.png)

```wl
Out[]= {0, 0, 0}
```

![0taqrc40bml7k](img/0taqrc40bml7k.png)

```wl
Out[]= 0
```

```wl
Out[]= 0
```

### Redundancies in EoS

```wl
In[]:= EE = ThermoConstraints[\[Delta]\[ScriptCapitalE], \[Delta]\[ScriptCapitalP], \[Delta]\[ScriptCapitalX]aa, \[Delta]\[ScriptCapitalX]\[Omega]\[Omega], \[Delta]\[ScriptCapitalX]a\[Omega], \[Delta]\[ScriptCapitalQ]\[ScriptL], \[Delta]\[ScriptCapitalN], \[Delta]\[ScriptCapitalN]\[ScriptL]] + O[\[ScriptE]]^4;
```

![1tu15k1ofddps](img/1tu15k1ofddps.png)

![0n2i3gwbuynxg](img/0n2i3gwbuynxg.png)

![0h18j5sfsed8d](img/0h18j5sfsed8d.png)

![11k6py565x7o7](img/11k6py565x7o7.png)

```wl
In[]:= \[Delta]EOS = EOSConstructor4 @@ \[Delta]Consti // RedundantSolser // Simplify
```

![19bxs7w5uvyqv](img/19bxs7w5uvyqv.png)

### Entropy current redundancy

![0c9i909mc57kx](img/0c9i909mc57kx.png)

![0na0akq06c1z9](img/0na0akq06c1z9.png)

![1li1zgnhcdv6u](img/1li1zgnhcdv6u.png)

### Invariants

![1jhmrn06pssh4](img/1jhmrn06pssh4.png)

```wl
In[]:= \[ScriptCapitalI]2[\[Delta]EOS]
 \[ScriptCapitalI]4[\[Delta]EOS]
```

```wl
Out[]= 0
```

```wl
Out[]= 0
```

```wl
In[]:= \[Delta]\[ScriptCapitalF] = (\[Delta]\[ScriptCapitalP] + vw2 T0^2 \[Delta]\[ScriptCapitalX]\[Omega]\[Omega]) /. ScaleVars // # + O[\[ScriptE]]^7 & // Normal // Simplify;
```

```wl
In[]:= \[ScriptCapitalI]2[\[Delta]\[ScriptCapitalF]]
 \[ScriptCapitalI]4[\[Delta]\[ScriptCapitalF]]
```

```wl
Out[]= 0
```

```wl
Out[]= 0
```

### Fixing redundancy for conformal EoS

![0rv303cdzfo10](img/0rv303cdzfo10.png)

![1sbe5mlx5t9fk](img/1sbe5mlx5t9fk.png)

```wl
In[]:= \[Delta]C = \[Delta]\[ScriptCapitalF]Conf + O[\[ScriptE]]^11 // Simplify;
```

![08nc5qf3co2bv](img/08nc5qf3co2bv.png)

```wl
In[]:= \[Delta]C /. sol2 /. sol4 /. sol6 /. sol8 /. sol10 // Simplify
```

![08oybofu52yp2](img/08oybofu52yp2.png)

### Thermodynamic Relations

![17zmxpr0tpuvo](img/17zmxpr0tpuvo.png)

![0xnbum8dhqfdd](img/0xnbum8dhqfdd.png)

#### Generic Expressions

![1tucahyu4lxb9](img/1tucahyu4lxb9.png)

```wl
In[]:= \[ScriptCapitalF]C = \[ScriptCapitalP]C + vw2 T0^2 \[ScriptCapitalX]\[Omega]\[Omega]C;
```

```wl
In[]:= GeneralConsti = {\[ScriptCapitalE]C, \[ScriptCapitalP]C, \[ScriptCapitalX]aaC, \[ScriptCapitalX]\[Omega]\[Omega]C, \[ScriptCapitalX]a\[Omega]C, \[ScriptCapitalQ]\[ScriptL]C, \[ScriptCapitalN]C, \[ScriptCapitalN]\[ScriptL]C};
```

```wl
In[]:= \[ScriptCapitalI]2[\[ScriptCapitalP]C + vw2 T0^2 \[ScriptCapitalX]\[Omega]\[Omega]C]
 \[ScriptCapitalI]4[\[ScriptCapitalP]C + vw2 T0^2 \[ScriptCapitalX]\[Omega]\[Omega]C]
```

![1qkv58sz4sa6s](img/1qkv58sz4sa6s.png)

![1snl37jkvbxzq](img/1snl37jkvbxzq.png)

#### Imposing Conservation Identities

![1c8ezgda1kq5o](img/1c8ezgda1kq5o.png)

```wl
Out[]= {0, 0}
```

#### Deriving Thermodynamic Relations

```wl
In[]:= ThermoEEC = ThermoConstraints @@ (GeneralConsti + \[Delta]Consti) /. thermoC // # + O[\[ScriptE]]^4 &;
```

![003cx0bxqk9xv](img/003cx0bxqk9xv.png)

```wl
In[]:= ClearAll[ThermodynamicRelations]
```

![1t6zryj2exs7y](img/1t6zryj2exs7y.png)

```wl
In[]:= ThermodynamicRelationsAns - (ThermodynamicRelations @@ GeneralConsti) // FullSimplify
```

```wl
Out[]= {0, 0, 0, 0, 0}
```

## Massless Dirac Fermions

### Expressions

![00kirumxtntj1](img/00kirumxtntj1.png)

```wl
In[]:= FermionConsti = {\[ScriptCapitalE]Fermion, \[ScriptCapitalP]Fermion, \[ScriptCapitalX]aaFermion, \[ScriptCapitalX]\[Omega]\[Omega]Fermion, \[ScriptCapitalX]a\[Omega]Fermion, \[ScriptCapitalQ]\[ScriptL]Fermion, \[ScriptCapitalN]Fermion, \[ScriptCapitalN]\[ScriptL]Fermion};
```

```wl
In[]:= ConservationIdentity @@ FermionConsti // Simplify
 TraceIdentity @@ FermionConsti // Simplify
```

```wl
Out[]= {0, 0}
```

```wl
Out[]= 0
```

```wl
In[]:= ThermodynamicRelations @@ FermionConsti // Simplify
```

```wl
Out[]= {0, 0, 0, 0, 0}
```

### Equation of state

```wl
In[]:= ThermoEEFermion = ThermoConstraints @@ (FermionConsti + \[Delta]Consti) // # + O[\[ScriptE]]^4 &;
```

![1ukc80ycd9t4f](img/1ukc80ycd9t4f.png)

![13enwsgoy3jrm](img/13enwsgoy3jrm.png)

```wl
In[]:= EOSFermion = EOSConstructor4 @@ (FermionConsti + \[Delta]Consti) // FermionSolser // # + O[\[ScriptE]]^5 &
```

![0gg3wxu9iob4q](img/0gg3wxu9iob4q.png)

```wl
In[]:= \[ScriptCapitalI]2[EOSFermion]
 \[ScriptCapitalI]4[EOSFermion]
```

![0200p7w9sl7k9](img/0200p7w9sl7k9.png)

![00ab0dgpa2m0o](img/00ab0dgpa2m0o.png)

## Massive Dirac Fermions

### Expressions

![0i01u0a2qkddb](img/0i01u0a2qkddb.png)

![1k0me1r498b1e](img/1k0me1r498b1e.png)

```wl
In[]:= MFermionConsti = {\[ScriptCapitalE]MFermion, \[ScriptCapitalP]MFermion, \[ScriptCapitalX]aaMFermion, \[ScriptCapitalX]\[Omega]\[Omega]MFermion, \[ScriptCapitalX]a\[Omega]MFermion, \[ScriptCapitalQ]\[ScriptL]MFermion, \[ScriptCapitalN]MFermion, \[ScriptCapitalN]\[ScriptL]MFermion};
```

```wl
In[]:= ConservationIdentity @@ MFermionConsti // FullSimplify
 TraceIdentity @@ MFermionConsti // FullSimplify
```

```wl
Out[]= {0, 0}
```

![1g8u6x3q29plj](img/1g8u6x3q29plj.png)

```wl
In[]:= ThermodynamicRelations @@ MFermionConsti // FullSimplify
```

```wl
Out[]= {0, 0, 0, 0, 0}
```

### Limit check (no need to evaluate)

```wl
In[]:= MasslessExpand[F_] := Sum[(-1)^(n + 1) F + O[m] // # + O[\[Nu]]^10 & // Normal, {n, 1, \[Infinity]}, Regularization -> "Dirichlet"]
```

```wl
In[]:= MasslessExpand[\[ScriptCapitalE]MFermion] - \[ScriptCapitalE]Fermion /. ScaleVars // # + O[\[ScriptE]]^4 & // Expand // Simplify
 MasslessExpand[\[ScriptCapitalP]MFermion] - \[ScriptCapitalP]Fermion /. ScaleVars // # + O[\[ScriptE]]^4 & // Expand // Simplify
 MasslessExpand[\[ScriptCapitalX]aaMFermion] - \[ScriptCapitalX]aaFermion /. ScaleVars // # + O[\[ScriptE]]^2 & // Expand // Simplify
 MasslessExpand[\[ScriptCapitalX]\[Omega]\[Omega]MFermion] - \[ScriptCapitalX]\[Omega]\[Omega]Fermion /. ScaleVars // # + O[\[ScriptE]]^2 & // Expand // Simplify
 MasslessExpand[\[ScriptCapitalX]a\[Omega]MFermion] - \[ScriptCapitalX]a\[Omega]Fermion /. ScaleVars // # + O[\[ScriptE]]^2 & // Expand // Simplify
 MasslessExpand[\[ScriptCapitalQ]\[ScriptL]MFermion] - \[ScriptCapitalQ]\[ScriptL]Fermion /. ScaleVars // # + O[\[ScriptE]]^2 & // Expand // Simplify
 MasslessExpand[\[ScriptCapitalN]MFermion] - \[ScriptCapitalN]Fermion /. ScaleVars // # + O[\[ScriptE]]^4 & // Expand // Simplify
 MasslessExpand[\[ScriptCapitalJ]\[ScriptL]MFermion] - \[ScriptCapitalJ]\[ScriptL]Fermion /. ScaleVars // # + O[\[ScriptE]]^2 & // Expand // Simplify
```

![1gqti664n23je](img/1gqti664n23je.png)

![1xcsdpkdsfm3g](img/1xcsdpkdsfm3g.png)

![0tjj3bi9nw3tz](img/0tjj3bi9nw3tz.png)

![01zjy68rwi7ug](img/01zjy68rwi7ug.png)

![0zwda4e6569sc](img/0zwda4e6569sc.png)

![0lelji6tkc4de](img/0lelji6tkc4de.png)

![11eh38r0spll6](img/11eh38r0spll6.png)

![0zbk9ij5centf](img/0zbk9ij5centf.png)

### Equation of state

```wl
In[]:= ThermoEEMFermion = ThermoConstraints @@ (MFermionConsti + \[Delta]Consti) // # + O[\[ScriptE]]^2 &;
```

![1l1scbsemrk2l](img/1l1scbsemrk2l.png)

|  |  |  |  |  |
| - | - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

```wl
In[]:= EOSMFermion = EOSConstructor2 @@ (MFermionConsti + \[Delta]Consti) // MFermionSolser // FullSimplify
```

![10fni5h6mfdf1](img/10fni5h6mfdf1.png)

```wl
In[]:= I2Fermion = EOSMFermion // \[ScriptCapitalI]2 // FullSimplify
```

![0kl4xrz9qq0u2](img/0kl4xrz9qq0u2.png)

## Massless Scalar Fields

### Expressions

![18589v4nyw8br](img/18589v4nyw8br.png)

```wl
In[]:= ScalarConsti = {\[ScriptCapitalE]Scalar, \[ScriptCapitalP]Scalar, \[ScriptCapitalX]aaScalar, \[ScriptCapitalX]\[Omega]\[Omega]Scalar, \[ScriptCapitalX]a\[Omega]Scalar, \[ScriptCapitalQ]\[ScriptL]Scalar, \[ScriptCapitalN]Scalar, \[ScriptCapitalN]\[ScriptL]Scalar};
```

```wl
In[]:= ConservationIdentity @@ ScalarConsti // Simplify
 TraceIdentity @@ ScalarConsti // FullSimplify
```

```wl
Out[]= {0, 0}
```

![19o3snrbhanvk](img/19o3snrbhanvk.png)

```wl
In[]:= ThermodynamicRelations @@ ScalarConsti
```

```wl
Out[]= {0, 0, 0, 0, 0}
```

### Equation of state

```wl
In[]:= ThermoEEScalar = ThermoConstraints @@ (ScalarConsti + \[Delta]Consti) // # + O[\[ScriptE]]^4 &;
```

![1gzevcnqwkaue](img/1gzevcnqwkaue.png)

![08mw97hwe91w7](img/08mw97hwe91w7.png)

```wl
In[]:= EOSScalar = EOSConstructor4 @@ (ScalarConsti + \[Delta]Consti) // ScalarSolser // # + O[\[ScriptE]]^5 &
```

![04p1052lhpref](img/04p1052lhpref.png)

```wl
In[]:= \[ScriptCapitalI]2[EOSScalar]
 \[ScriptCapitalI]4[EOSScalar]
```

![1u3yohpm6v7sb](img/1u3yohpm6v7sb.png)

![1af4xkjh0wekj](img/1af4xkjh0wekj.png)

## Massive Scalar Fields

### Scalar Expressions

![1hq5j67nc9z3h](img/1hq5j67nc9z3h.png)

![1otco0ekq0rig](img/1otco0ekq0rig.png)

```wl
In[]:= MScalarConsti = {\[ScriptCapitalE]MScalar, \[ScriptCapitalP]MScalar, \[ScriptCapitalX]aaMScalar, \[ScriptCapitalX]\[Omega]\[Omega]MScalar, \[ScriptCapitalX]a\[Omega]MScalar, \[ScriptCapitalQ]\[ScriptL]MScalar, \[ScriptCapitalN]MScalar, \[ScriptCapitalN]\[ScriptL]MScalar};
```

```wl
In[]:= ConservationIdentity @@ MScalarConsti // FullSimplify
 TraceIdentity @@ MScalarConsti // FullSimplify
```

```wl
Out[]= {0, 0}
```

![1el1xcg4o9fe7](img/1el1xcg4o9fe7.png)

```wl
In[]:= ThermodynamicRelations @@ MScalarConsti // FullSimplify
```

```wl
Out[]= {0, 0, 0, 0, 0}
```

### Limit check (no need to evaluate)

```wl
In[]:= MasslessExpandBoson[F_] := Sum[F + O[m] // # + O[\[Nu]]^10 & // Normal, {n, 1, \[Infinity]}, Regularization -> "Dirichlet"]
```

```wl
In[]:= MasslessExpandBoson[\[ScriptCapitalE]MScalar] - \[ScriptCapitalE]Scalar + O[\[ScriptE]]^4 /. ScaleVars // Expand // FullSimplify
 MasslessExpandBoson[\[ScriptCapitalP]MScalar] - \[ScriptCapitalP]Scalar + O[\[ScriptE]]^4 /. ScaleVars // Expand // FullSimplify
 MasslessExpandBoson[\[ScriptCapitalX]aaMScalar] - \[ScriptCapitalX]aaScalar /. ScaleVars // # + O[\[ScriptE]]^2 & // Expand // Simplify
 MasslessExpandBoson[\[ScriptCapitalX]\[Omega]\[Omega]MScalar] - \[ScriptCapitalX]\[Omega]\[Omega]Scalar /. ScaleVars // # + O[\[ScriptE]]^2 & // Expand // Simplify
 MasslessExpandBoson[\[ScriptCapitalX]a\[Omega]MScalar] - \[ScriptCapitalX]a\[Omega]Scalar /. ScaleVars // # + O[\[ScriptE]]^2 & // Expand // Simplify
 MasslessExpandBoson[\[ScriptCapitalQ]\[ScriptL]MScalar] - \[ScriptCapitalQ]\[ScriptL]Scalar /. ScaleVars // # + O[\[ScriptE]]^2 & // Expand // Simplify
 MasslessExpandBoson[\[ScriptCapitalN]MScalar] - \[ScriptCapitalN]Scalar /. ScaleVars // # + O[\[ScriptE]]^6 & // Expand // Simplify
 MasslessExpandBoson[\[ScriptCapitalN]\[ScriptL]MScalar] - \[ScriptCapitalN]\[ScriptL]Scalar /. ScaleVars // # + O[\[ScriptE]]^6 & // Expand // Simplify
```

![1xbgc05bla63t](img/1xbgc05bla63t.png)

![192783wome3ln](img/192783wome3ln.png)

![0umlhq3bc8aqh](img/0umlhq3bc8aqh.png)

![17lo4ux84ynaj](img/17lo4ux84ynaj.png)

![10wg9ip71n2oj](img/10wg9ip71n2oj.png)

![1cggogsle0dr3](img/1cggogsle0dr3.png)

![1xbcpxgfi5q0b](img/1xbcpxgfi5q0b.png)

![1wzq4d9zrg5cl](img/1wzq4d9zrg5cl.png)

![1mx1baqi8v3pr](img/1mx1baqi8v3pr.png)

![1efhqkgpv695k](img/1efhqkgpv695k.png)

![0othkcwhajtle](img/0othkcwhajtle.png)

![11mfaeaf85d43](img/11mfaeaf85d43.png)

![038b5j687l6ms](img/038b5j687l6ms.png)

![1k81yf65lvqtz](img/1k81yf65lvqtz.png)

![0s55ey0l08i2d](img/0s55ey0l08i2d.png)

![1usaiw5aaiupe](img/1usaiw5aaiupe.png)

![0omuxa3n9578q](img/0omuxa3n9578q.png)

![01goe87hjdoko](img/01goe87hjdoko.png)

![0s6f5jfr4v3i5](img/0s6f5jfr4v3i5.png)

```wl
Out[]= 0
```

![0dsq1vdogxvwt](img/0dsq1vdogxvwt.png)

### Equation of state

```wl
In[]:= ThermoEEMScalar = ThermoConstraints @@ (MScalarConsti + \[Delta]Consti) // # + O[\[ScriptE]]^2 &;
```

![10hpm1plpvmpp](img/10hpm1plpvmpp.png)

|  |  |  |  |  |
| - | - | - | - | - |
| -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- | -SeriesData- |

```wl
In[]:= EOSMScalarVal = EOSConstructor2 @@ (MScalarConsti + \[Delta]Consti) // MScalarSolser // # + O[\[ScriptE]]^3 & // FullSimplify
```

![0ppr51u5tlxew](img/0ppr51u5tlxew.png)

![12ejq2ypcjoza](img/12ejq2ypcjoza.png)

![0vjcft9z5ria5](img/0vjcft9z5ria5.png)

```wl
In[]:= \[ScriptCapitalI]2[EOSMScalar] /. \[ScriptE] -> 1 // FullSimplify
```

![0h37sbnsovyww](img/0h37sbnsovyww.png)

## Export

```wl
In[]:= Export["notebook.md", EvaluationNotebook[]]
```