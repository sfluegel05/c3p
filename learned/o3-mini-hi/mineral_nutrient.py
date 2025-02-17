"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: Mineral nutrient
Definition: A mineral that is an inorganic nutrient which must be ingested and absorbed 
in adequate amounts to satisfy a wide range of essential metabolic and/or structural functions
in the human body.
We try to identify salts that contain a metal ion (or an inorganic complex) from a pre‐defined set,
and reject structures where the counterion fragments are “organic” (e.g. contain aromatic rings or unsaturated C=C bonds)
or where the metal ion is not “nutrient‐like” (for example if the metal has a negative formal charge or is coordinated by –OH groups).
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    Our heuristic uses three ideas:
      (1) At least one fragment must contain a metal from a defined set.
      (2) If the metal appears as a free ion (i.e. a single-atom fragment) then it must show a positive formal charge.
      (3) If the metal is in a multicomponent (complex) fragment, then that fragment must be inorganic 
          (i.e. it should not contain any carbon atoms) OR if it does, the ligand should be very simple
          (e.g. a saturated carboxylate, not an unsaturated one) and the metal’s coordination environment should not include –OH groups.
      (4) In addition, if any fragment (that is not just water) contains carbon and shows either aromatic rings or non-carbonyl C=C bonds,
          we assume it is part of an organic dye or complex salt and we reject the candidate.
      (5) In a special case we disallow e.g. aluminium chloride type complexes.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a mineral nutrient, False otherwise.
        str: Reason for the classification.
        
    Note:
        If the task is too hard to decide the classification, the function may return (None, None).
    """
    # Try to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # It helps to have explicit hydrogens to check for –OH groups.
    mol = Chem.AddHs(mol)
    
    # Define a set of metal elements that we consider candidates.
    metal_set = {"Pd", "K", "Fe", "Cs", "Ca", "Zn", "Al", "Mg", "Sb", "Ba", "Na", "La"}
    
    # Helper function: check if an atom (of O) is –OH: i.e. it has at least one neighboring hydrogen.
    def is_hydroxyl(oxygen):
        for nbr in oxygen.GetNeighbors():
            if nbr.GetSymbol() == "H":
                return True
        return False
    
    # Helper function: for a given fragment (mol), check if it is "organic" in a disqualifying way.
    def fragment_is_organic(frag):
        # If the fragment contains carbon, then check for aromatic rings.
        if any(atom.GetSymbol() == "C" for atom in frag.GetAtoms()):
            # If any atom is aromatic, disqualify.
            if any(atom.GetIsAromatic() for atom in frag.GetAtoms()):
                return True
            # Also, check each bond: if there’s any double bond between two carbons,
            # and it is not part of a carbonyl (i.e. C=O), then disqualify.
            for bond in frag.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    a1 = bond.GetBeginAtom()
                    a2 = bond.GetEndAtom()
                    if a1.GetSymbol() == "C" and a2.GetSymbol() == "C":
                        # For a carbonyl, one of these carbons should be double-bonded to oxygen.
                        is_carbonyl = False
                        for nb in a1.GetNeighbors():
                            if nb.GetSymbol() == "O":
                                for b in frag.GetBonds():
                                    if (b.GetBeginAtom() == a1 and b.GetEndAtom() == nb) or (b.GetBeginAtom() == nb and b.GetEndAtom() == a1):
                                        if b.GetBondType() == Chem.BondType.DOUBLE:
                                            is_carbonyl = True
                        for nb in a2.GetNeighbors():
                            if nb.GetSymbol() == "O":
                                for b in frag.GetBonds():
                                    if (b.GetBeginAtom() == a2 and b.GetEndAtom() == nb) or (b.GetBeginAtom() == nb and b.GetEndAtom() == a2):
                                        if b.GetBondType() == Chem.BondType.DOUBLE:
                                            is_carbonyl = True
                        if not is_carbonyl:
                            return True
        return False

    # Get fragments (using the dot notation in SMILES gives separate fragments)
    frags = Chem.rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    
    metal_found = False
    reasons = []
    disqualifying_fragment = None

    # Check each fragment.
    for frag in frags:
        # Get the set of atomic symbols present in this fragment.
        atom_syms = [atom.GetSymbol() for atom in frag.GetAtoms()]
        
        # If fragment is water (H2O) we ignore it.
        if len(frag.GetAtoms()) == 3 and sorted(atom_syms) == ["H", "H", "O"]:
            continue
        
        # Check if the fragment appears organic (has carbon and disqualifying features)
        if fragment_is_organic(frag):
            disqualifying_fragment = Chem.MolToSmiles(frag)
            reasons.append("Contains an organic fragment (aromatic or unsaturated C=C): " + disqualifying_fragment)
        
        # Now examine atoms in this fragment to see if any is a metal from our set.
        # We require that either the metal appears as a free ion (a one-atom fragment with positive charge)
        # or it is embedded in an inorganic fragment.
        for atom in frag.GetAtoms():
            sym = atom.GetSymbol()
            if sym not in metal_set:
                continue

            # Rule A: if the fragment is a single atom, we require a positive formal charge.
            if frag.GetNumAtoms() == 1:
                if atom.GetFormalCharge() > 0:
                    metal_found = True
                    reasons.append(f"Found isolated metal ion: {sym}{atom.GetFormalCharge()}")
            else:
                # For metal atoms in a fragment with more than one atom, first reject if the formal charge is negative.
                if atom.GetFormalCharge() < 0:
                    continue
                # If there is at least one carbon in the fragment, we treat it as organic/inorganic hybrid.
                has_carbon = any(a.GetSymbol() == "C" for a in frag.GetAtoms())
                if not has_carbon:
                    # For pure inorganic fragments (no carbon) we also check that no –OH groups (oxygen with H)
                    # are directly bonded to this metal. (This helps to remove cases like antimonic acid.)
                    rejects_due_to_OH = False
                    for nbr in atom.GetNeighbors():
                        if nbr.GetSymbol() == "O" and is_hydroxyl(nbr):
                            rejects_due_to_OH = True
                            break
                    # Also, as a special case disallow aluminium complexes where all the ligands are chlorides.
                    if sym == "Al":
                        all_cl = True
                        for nbr in atom.GetNeighbors():
                            if nbr.GetSymbol() != "Cl":
                                all_cl = False
                                break
                        if all_cl:
                            rejects_due_to_OH = True  # use this flag to disqualify
                    if not rejects_due_to_OH:
                        metal_found = True
                        reasons.append(f"Found inorganic metal in non-organic fragment: {sym}{'' if atom.GetFormalCharge()==0 else atom.GetFormalCharge()}")
                else:
                    # In a fragment that contains carbon, if we find the metal we only accept it
                    # if the organic part is very simple. For our heuristic, we allow a carbon-containing counterion
                    # if it is a small carboxylate (i.e. if it does NOT contain any C=C bonds outside of C=O)
                    # We already flagged fragments with aromatic rings or unsaturated C=C bonds as disqualifying.
                    if not fragment_is_organic(frag):
                        # For example, in magnesium dipropionate the metal [Mg++] is found in a fragment separated as a free ion;
                        # the fatty acid fragments (if they are saturated) are acceptable counterions.
                        metal_found = True
                        reasons.append(f"Found metal in a carbon-containing but simple fragment: {sym}{'' if atom.GetFormalCharge()==0 else atom.GetFormalCharge()}")

    # If any disqualifying (organic) fragment was found, then we do not want to classify this as a mineral nutrient.
    if disqualifying_fragment:
        return False, "Disqualifying organic fragment detected: " + disqualifying_fragment

    if metal_found:
        reason_str = " ; ".join(reasons)
        return True, "Found metal nutrient ion(s): " + reason_str
    else:
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
        # False positive examples:
        "tetrachloroaluminate(1-)": "Cl[Al-](Cl)(Cl)Cl",
        "trypan blue": "[Na+].[Na+].[Na+].[Na+].Cc1cc(ccc1\\N=N\\c1c(O)c2c(N)cc(cc2cc1S([O-])(=O)=O)S([O-])(=O)=O)-c1ccc(\\N=N\\c2c(O)c3c(N)cc(cc3cc2S([O-])(=O)=O)S([O-])(=O)=O)c(C)c1",
        "NIR-3 dye": "[K+].[K+].[H]C(=CC([H])=CC([H])=C1N(CCCCS([O-])(=O)=O)c2ccc(cc2C1(C)C)C(O)=O)C=C([H])C1=[N+](CCCCS([O-])(=O)=O)c2ccc(cc2C1(C)C)S([O-])(=O)=O",
        "antimonic acid": "[H]O[Sb](=O)(O[H])O[H]",
        "prednisolone sodium succinate": "[Na+].[H][C@@]12CC[C@](O)(C(=O)COC(=O)CCC([O-])=O)[C@@]1(C)C[C@H](O)[C@@]1([H])[C@@]2([H])CCC2=CC(=O)C=C[C@]12C",
        "zincide": "[Zn-]",
        "aluminium trichloride hexahydrate": "[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].Cl[Al](Cl)Cl",
        "Cefditoren Sodium Salt": "S1C2N(C(C(=O)[O-])=C(C1)C=CC=3SC=NC3C)C(=O)C2NC(=O)C(=NOC)C=4N=C(N)SC4.[Na+]",
        "hexaaquamagnesium(2+)": "[H][O]([H])[Mg++]([O]([H])[H])([O]([H])[H])([O]([H])[H])([O]([H])[H])[O]([H])[H]",
        "acid red 29": "[Na+].[Na+].Oc1cc(cc2cc(c(N=Nc3ccccc3)c(O)c12)S([O-])(=O)=O)S([O-])(=O)=O",
        "aluminium hydroxide": "[H]O[Al](O[H])O[H]",
        "3-vinylbacteriochlorophyllide a(1-)": "CC[C@@H]1[C@@H](C)C2=Cc3c(C=C)c(C)c4C=C5[C@@H](C)[C@H](CCC([O-])=O)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31",
        "DY-634": "[Na+].[Na+].[Na+].[H]C(=C([H])c1cc([o+]c2cc(ccc12)N(CCCS([O-])(=O)=O)CCCS([O-])(=O)=O)C(C)(C)C)C([H])=C1N(CCCC(O)=O)c2ccc(cc2C1(C)CCCS([O-])(=O)=O)S([O-])(=O)=O",
        "ferrienterobactin": "[H][C@]12COC(=O)[C@]3([H])COC(=O)[C@]([H])(COC1=O)[NH2+]C(=O)c1cccc4O[Fe-3]56(Oc7cccc(C(=O)[NH2+]2)c7O5)(Oc2cccc(C(=O)[NH2+]3)c2O6)Oc14",
        "phenoxymethylpenicillin potassium": "[K+].[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)COc1ccccc1)C([O-])=O",
        "Aceglutamide aluminum": "NC(CC[C@@H](C(=O)[O-])NC(C)=O)=O.[Al+3].[Al+3].[Al+3].[OH-].[OH-].[OH-].[OH-].NC(CC[C@@H](C(=O)[O-])NC(C)=O)=O.NC(CC[C@@H](C(=O)[O-])NC(C)=O)=O.NC(CC[C@@H](C(=O)[O-])NC(C)=O)=O.NC(CC[C@@H](C(=O)[O-])NC(C)=O)=O",
        "sodium hexacyanoferrate(4-)": "[Na+].[Na+].[Na+].[Na+].N#C[Fe-4](C#N)(C#N)(C#N)(C#N)C#N",
        "7(1)-hydroxychlorophyll a(1-)": "CCC1=C(CO)C2=[N+]3C(=CC4=[N+]5C(=CC=6N7C=8C(=C9[N+](=C(C1)[C@H]([C@@H]9CCC(O)=O)C)[Mg-2]753)CC(C8C6CC)=O)C(=C4C)CC)C(=C2C)C=C",
        "Calcium linoleate": "C(CC([O-])=O)CCCCC/C=C\\C/C=C\\CCCCC.[Ca+2].C(CC([O-])=O)CCCCC/C=C\\C/C=C\\CCCCC"
    }
    
    for name, smi in test_smiles.items():
        result, reason = is_mineral_nutrient(smi)
        print(f"{name}: {result}, {reason}")