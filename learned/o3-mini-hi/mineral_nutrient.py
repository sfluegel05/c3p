"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: Mineral nutrient
Definition: A mineral that is an inorganic nutrient which must be ingested and
absorbed in adequate amounts to satisfy a wide range of essential metabolic and/or 
structural functions in the human body.
Heuristic idea:
  - At least one fragment must supply an acceptable metal ion (from a defined set)
    either as a free ion (one‐atom fragment with positive charge) or as part of a clearly inorganic fragment.
  - Additionally, if a free metal ion is detected then every other fragment (counterion)
    must not look “organic” in a way that suggests it is a small, water‐soluble organic acid
    (for instance sodium acetate is disallowed even though [Na+] is isolated, while caesium formate is allowed).
  - We flag fragments as “organic” if they have carbon and either (a) contain an aromatic atom 
    or (b) show a C=C bond that is not part of a carbonyl.
  
Because the boundaries are not sharp this implementation is only one possible heuristic.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    
    Our heuristic has three main parts:
      (1) The molecule must supply at least one acceptable metal ion.
          Acceptable metal ions are those whose elemental symbol
          appears in our allowed metal set and which (a) appear as an isolated ion 
          (i.e. the only atom in that fragment with a positive formal charge) or (b)
          appear in an entirely inorganic fragment (i.e. no carbon atoms).
      (2) Any fragment that does not contain a metal (the counterion) should be checked.
          Here we allow fragments with no carbon or – if they do contain carbon – only those
          that are either very simple (e.g. a formate, with one carbon, or a propionate, with three carbons)
          or very lipophilic (in our crude test: having at least 8 carbons).
      (3) Independently, if any fragment (other than water) is overall “organic”
          (for example having aromatic rings or non‐carbonyl C=C bonds), we flag the entire molecule as disqualified.
    
    If these criteria are met then we return True plus a reason string; otherwise False with a reason.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple of classification and explanation.
    """
    # Parse SMILES and add explicit hydrogens (to help assess –OH groups later)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Define allowed metal elements
    metal_set = {"Pd", "K", "Fe", "Cs", "Ca", "Zn", "Al", "Mg", "Sb", "Ba", "Na", "La"}
    
    # ---------------------------------------------------------------------
    # Helper functions
    # ---------------------------------------------------------------------
    def is_water(frag):
        # Water should have exactly two H and one O.
        at_syms = sorted([atom.GetSymbol() for atom in frag.GetAtoms()])
        return len(frag.GetAtoms()) == 3 and at_syms == ["H", "H", "O"]
    
    def fragment_is_organic(frag):
        """
        Flag a fragment as organic if it contains a carbon and either any aromatic atoms,
        or any double bond between carbons that is not part of a carbonyl.
        """
        if any(atom.GetSymbol() == "C" for atom in frag.GetAtoms()):
            if any(atom.GetIsAromatic() for atom in frag.GetAtoms()):
                return True
            for bond in frag.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    a1 = bond.GetBeginAtom()
                    a2 = bond.GetEndAtom()
                    if a1.GetSymbol()=="C" and a2.GetSymbol()=="C":
                        # Check if one of these is in a carbonyl group.
                        def is_carbonyl(atom):
                            for nbr in atom.GetNeighbors():
                                b = frag.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                                if nbr.GetSymbol() == "O" and b and b.GetBondType() == Chem.BondType.DOUBLE:
                                    return True
                            return False
                        if not (is_carbonyl(a1) or is_carbonyl(a2)):
                            return True
        return False
    
    def carbon_count(frag):
        return sum(1 for atom in frag.GetAtoms() if atom.GetSymbol() == "C")
    
    # For a non‐metal fragment (a counterion) we allow it if either:
    #  — It contains no carbon, or
    #  — It does contain carbon but appears “simple”: either it is very small (1 or 3 C: formate or propionate)
    #     or—if larger (8 or more carbons)—we assume it is a fatty‐acid type ligand.
    def allowed_counterion(frag):
        cnt = carbon_count(frag)
        if cnt == 0:
            return True
        if cnt in [1, 3]:
            return True
        if cnt >= 8:
            # But if it is obviously organic (for example showing aromatic atoms or non‐carbonyl double bonds),
            # then reject.
            if fragment_is_organic(frag):
                return False
            return True
        # Otherwise (2, 4–7 carbons) we do not allow it.
        return False
    
    # ---------------------------------------------------------------------
    # Main analysis
    # ---------------------------------------------------------------------
    # Split the molecule into fragments (based on dot notation)
    frags = Chem.rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    
    free_metal_found = False
    inorganic_metal_found = False
    metal_reasons = []
    counterions_ok = True
    counterion_reasons = []
    overall_disqual = False
    dq_reason = ""
    
    # (A) First pass: if any fragment (other than water) appears organic in a disqualifying way, flag molecule.
    for frag in frags:
        if is_water(frag):
            continue
        # We check fragments that do not contain a metal from our allowed set.
        if not any(atom.GetSymbol() in metal_set for atom in frag.GetAtoms()):
            if any(atom.GetSymbol()=="C" for atom in frag.GetAtoms()):
                if fragment_is_organic(frag):
                    overall_disqual = True
                    dq_reason = "Disqualifying organic fragment detected: " + Chem.MolToSmiles(frag)
                    break
    
    # (B) Now look for metal‐containing fragments.
    for frag in frags:
        if is_water(frag):
            continue
        for atom in frag.GetAtoms():
            sym = atom.GetSymbol()
            if sym not in metal_set:
                continue
            # If the fragment consists only of the metal atom, then we require a positive formal charge.
            if frag.GetNumAtoms() == 1:
                if atom.GetFormalCharge() > 0:
                    free_metal_found = True
                    metal_reasons.append(f"Found isolated metal ion: {sym}{atom.GetFormalCharge()}")
            else:
                # In a multi‐atom fragment:
                # if no carbon is present, then we accept it as an inorganic complex.
                if carbon_count(frag) == 0:
                    if atom.GetFormalCharge() >= 0:
                        inorganic_metal_found = True
                        metal_reasons.append(f"Found inorganic metal in non-organic fragment: {sym}{'' if atom.GetFormalCharge()==0 else atom.GetFormalCharge()}")
                # If metal is in a fragment that contains carbon we consider it “complex”
                # and (by our simple heuristic) we reject that metal.
    
    # (C) Finally, if a free metal ion has been found, examine every other fragment (that does NOT contain a metal)
    # as a “counterion”. If any such fragment looks unacceptable (i.e. its structure is not allowed by allowed_counterion),
    # then we disqualify the whole molecule.
    for frag in frags:
        if is_water(frag):
            continue
        # Skip fragments that already contain a metal.
        if any(atom.GetSymbol() in metal_set for atom in frag.GetAtoms()):
            continue
        if any(atom.GetSymbol()=="C" for atom in frag.GetAtoms()):
            if not allowed_counterion(frag):
                counterions_ok = False
                counterion_reasons.append("Disallowed counterion fragment: " + Chem.MolToSmiles(frag))
    
    # (D) Now summarize:
    if overall_disqual:
        return False, dq_reason
    # We require at least one acceptable metal ion.
    if free_metal_found or inorganic_metal_found:
        # In case the acceptable metal was present only as a free ion, then we also demand that no disallowed counterions are present.
        if free_metal_found and not counterions_ok:
            reason = ("Found free metal ion(s): " + " ; ".join(metal_reasons) +
                      " but disallowed counterion fragment(s): " + " ; ".join(counterion_reasons))
            return False, reason
        reason = "Found metal nutrient ion(s): " + " ; ".join(metal_reasons)
        return True, reason
    return False, "No acceptable metal nutrient ion found"


# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = {
      "Potassium hexachloropalladate(IV)": "[Pd-2](Cl)(Cl)(Cl)(Cl)(Cl)Cl.[K+].[K+]",
      "iron(3+) phosphate": "[Fe+3].[O-]P([O-])(=O)[O-]",
      "caesium formate": "[Cs+].[H]C([O-])=O",
      "calcium silicate": "[Ca++].[Ca++].[O-][Si]([O-])([O-])[O-]",
      "zinc nitrate": "[Zn++].[O-][N+]([O-])=O.[O-][N+]([O-])=O",
      "aluminium sulfate octadecahydrate": "O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.[Al+3].[Al+3].[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O",
      "magnesium distearate": "[Mg++].CCCCCCCCCCCCCCCCCC([O-])=O.CCCCCCCCCCCCCCCCCC([O-])=O",
      "antimony pentafluoride": "[Sb](F)(F)(F)(F)F",
      "caesium chloride": "[Cl-].[Cs+]",
      "Calcium hypochlorite": "Cl[O-].[Ca+2].Cl[O-]",
      "Lanthanum trichloride": "Cl[La](Cl)Cl",
      "magnesium phosphate": "P([O-])([O-])([O-])=O.P([O-])([O-])([O-])=O.[Mg+2].[Mg+2].[Mg+2]",
      "calcium dichloride": "[Cl-].[Cl-].[Ca++]",
      "trisodium phosphate": "[Na+].[Na+].[Na+].[O-]P([O-])([O-])=O",
      "barium sulfate": "[O-]S([O-])(=O)=O.[Ba+2]",
      "aluminium sulfate (anhydrous)": "[Al+3].[Al+3].[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O",
      "barium nitrate": "[Ba++].[O-][N+]([O-])=O.[O-][N+]([O-])=O",
      "calcium monohydroxide": "O[Ca]",
      "calcium sulfate": "[Ca++].[O-]S([O-])(=O)=O",
      "magnesium dipropionate": "[Mg++].CCC([O-])=O.CCC([O-])=O",
      "tricalcium bis(phosphate)": "[Ca++].[Ca++].[Ca++].[O-]P([O-])([O-])=O.[O-]P([O-])([O-])=O",
      "disodium hydrogenphosphate": "[Na+].[Na+].OP([O-])([O-])=O",
      "potassium chloride": "[Cl-].[K+]",
      "calcium dihydroxide": "[OH-].[OH-].[Ca++]",
      "magnesium dichloride": "[Mg++].[Cl-].[Cl-]",
      "magnesium sulfate": "[Mg++].[O-]S([O-])(=O)=O",
      "barium acetate": "[Ba++].CC([O-])=O.CC([O-])=O",
      "barium carbonate": "[Ba++].[O-]C([O-])=O",
      "calcium carbonate": "[Ca+2].C(=O)([O-])[O-]",
      "calcium difluoride": "[F-].[F-].[Ca++]",
      "potassium sulfate": "[K+].[K+].[O-]S([O-])(=O)=O",
      "calcium hydrogenphosphate": "[Ca++].[H]OP([O-])([O-])=O",
      # some false positives (expected to be rejected)
      "Aceglutamide aluminum": "NC(CC[C@@H](C(=O)[O-])NC(C)=O)=O.[Al+3].[Al+3].[Al+3].[OH-].[OH-].[OH-].[OH-].NC(CC[C@@H](C(=O)[O-])NC(C)=O)=O.NC(CC[C@@H](C(=O)[O-])NC(C)=O)=O.NC(CC[C@@H](C(=O)[O-])NC(C)=O)=O.NC(CC[C@@H](C(=O)[O-])NC(C)=O)=O",
      "sodium hexacyanoferrate(4-)": "[Na+].[Na+].[Na+].[Na+].N#C[Fe-4](C#N)(C#N)(C#N)(C#N)C#N",
      "sodium acetate": "[Na+].CC(=O)[O-]",
    }
    
    for name, smi in test_smiles.items():
        result, reason = is_mineral_nutrient(smi)
        print(f"{name}: {result}, {reason}")