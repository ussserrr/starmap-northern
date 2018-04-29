#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from collections import OrderedDict

from util import download_stars, extract_forms

# greek_alphabet = "αβγδεζηθικλμνξοπρστυφχψω"  # just for copying letters


# entry point (input data): dict with
#     keys: constellations' 3-letter code and
#     values: their forms as string of greek letters (stars of this constellation)
#             in concrete order to render each constellation
# sort dict by names of constellations to be sure that two tables corresponds to each other
constellations = OrderedDict(sorted(
    { 'UMa': ['η','ζ','ε','δ','γ','β','α','δ'],
      'UMi': ['ζ','β','γ','η','ζ','ε','δ','α'],
      'Her': ['ι','τ','φ','χ','σ','ι','θ','ρ','π','η','ζ','ε','β','γ','β','ε','μ','δ','α','δ','μ','ο','ξ','π','ξ','ο','109','110','113'],
      'Cas': ['ε','δ','γ','α','β'],
      'Dra': ['λ','κ','α','ι','θ','η','ζ','χ','τ','ε','δ','ξ','ν','β','γ','ξ'],
      'And': ['α','δ','β','μ','ν','μ','β','γ1'],
      'Gem': ['η','μ','ε','ν','ε','τ','θ','τ','α','τ','ι','υ','β','υ','κ','υ','δ','λ','ξ','λ','δ','ζ','γ'],
      'Lib': ['θ','γ','β','α','σ','γ'],  # Весы
      'Vir': ['μ','ι','κ','α','ζ','τ','109','τ','ζ','δ','ε','δ','γ','α','γ','η','ν'],
      'Ser': ['μ','ε','α','δ','β','κ','γ','β'],
      'Boo': ['θ','λ','γ','β','δ','μ','δ','ε','α','ζ','α','η','υ','η','α','ρ','γ'],
      'Leo': ['ε','μ','ζ','γ','η','α','θ','β','δ','γ'],
      'Cnc': ['χ','γ','ι','γ','δ','α','δ','β'],
      'Ori': ['κ','ζ','ε','δ','β','δ','γ','π3','π2','π1','π2','π3','π4','π5','π6','π5','π4','π3','γ','λ','α','ζ','α','μ','ξ','μ','ν','χ1','χ2','ξ'],
      'Tri': ['α','β','γ','α'],
      'Per': ['ζ','ξ','ε','δ','α','γ','η','γ','α','β','ρ'],
      'Cep': ['ι','β','α','ζ','ι','γ','β'],  # Цефей
      'Lyr': ['α','ζ','δ','γ','β','ζ'],
      'Cyg': ['β','η','γ','α','γ','δ','ι','κ','ι','δ','γ','ε','ζ'],
      'Psc': ['ι','λ','κ','γ','θ','ι','ome','d','δ','ε','ζ','μ','ν','ξ','α','ο','π','η','φ','σ','τ','υ','φ'],
      'Aql': ['γ','α','β','α','δ','λ','δ','η','θ','η','δ','ζ','ε'],
      'Lep': ['θ','η','ζ','α','β','γ','δ','α','μ','λ','ν','λ','μ','κ','ι','κ','μ','ε','β'],
      'CMa': ['ι','θ','γ','ι','α','ν2','β','ν2','ξ2','ν2','ο1','σ','ε','κ','ε','ζ','ε','σ','δ','ome','η','ome','δ','ο2','α'],
      'Crt': ['γ','β','α','δ','ε','θ','η','ζ','γ'],
      'Oph': ['c','η','ε','κ','α','β','η'],
      'Peg': ['ε','θ','ζ','ξ','α','β','η','π','η','β','μ','λ','ι','κ','ι','λ','μ','β','δ','γ','α'],
      'Aur': ['ι','ζ','α','β','θ','ι'],  # Возничий
      'Lyn': ['α','38','31','21','15','10','2'],
      'Aqr': ['ψ2','χ','φ','λ','η','ζ1','γ','α','θ','ι','θ','σ','τ','δ','τ','σ','θ','α','β','μ','ε'],
      'Sex': ['γ','α','β','ε','γ'],
      'Hya': ['ρ','η','σ','δ','ε','ρ','ζ','ome','θ','ι','τ1','α','υ1','υ2','λ','μ','ν','ξ','β','ψ','γ'],
      'Ant': ['α','η'],
      'Pyx': ['β','α','γ'],
      'Sct': ['α','β','ε','δ','γ','ζ','α'],
      'Cap': ['α','β','ψ','ome','θ','β','θ','ζ','ι','θ','ι','γ','δ'],
      'Tau': ['ο','λ','γ','θ1','α','ζ','α','ε','τ','β','τ','ε','δ3','δ2','δ3','η','δ3','γ'],
      'Cet': ['ξ1','ξ2','μ','λ','α','γ','ν','ξ2','ν','γ','δ','ο','ε','ρ','ζ','θ','η','β','ι','β','τ','σ','π','ε'],
#      'Vel': ['γ','λ','ψ','q','p','μ','φ','N','κ','δ','ο','γ'],
#      'Pup': ['τ','ν','π','ξ','ρ','ζ','σ'],
    }.items()))


if __name__ == '__main__':

    # download constellations
    stars_table = download_stars(constellations)
    stars_table.write('stars.html', format='ascii.html', overwrite=True)

    # create constellations forms
    constellations_forms = extract_forms(constellations, stars_table)
    constellations_forms.write('.constellations.html', format='ascii.html', overwrite=True)
