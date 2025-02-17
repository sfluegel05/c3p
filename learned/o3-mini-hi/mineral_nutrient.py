"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: Mineral nutrient
Definition: A mineral nutrient is defined as a mineral (from a set of allowed metals)
that is an inorganic nutrient essential for metabolic/structural function.
A heuristic is used:
  (A) There must be at least one acceptable metal‐containing fragment.
      Acceptable metal fragments are either isolated metal ions (one atom, positive charge) 
      or multi‐atom inorganic fragments that do not contain any carbon. 
      However, if a metal is coordinated to water (i.e. the only other atoms are O and H)
      then we treat it as a metal–aqua complex and disqualify that molecule.
  (B) Any fragment that does not contain an allowed metal is taken as a counterion.
      If it contains carbon then it must be “simple” (if a carboxylate) or—in the case of
      larger fragments (≥8 carbons)—must show no additional “organic” features (e.g. aromatic or non‐carbonyl double bonds).
Because the boundaries are not sharp, this is only one possible heuristic.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    
    Args:
       smiles (str): SMILES string.
       
    Returns:
       (bool, str): Tuple with overall classification and explanation.
    """
    # Parse the SMILES string. If it fails, return an error.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # For fragment evaluation, work with fragments as separate molecules.
    # (Don’t add explicit H for the purpose of fragment splitting.)
    frags = Chem.rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    
    # Define allowed metal symbols (as strings) and which ones are considered alkali.
    allowed_metals = {"Pd", "K", "Fe", "Cs", "Ca", "Zn", "Al", "Mg", "Sb", "Ba", "Na", "La"}
    alkali_metals = {"Na", "K", "Cs"}
    
    # ---------------- Helper functions ----------------
    
    def is_water(fragment):
        # Use molecular formula to detect water.
        # Note that water in SMILES may be "O" with implicit hydrogens.
        return Chem.CalcMolFormula(fragment) == "H2O"
    
    def removeHs(fragment):
        # Remove explicit hydrogens for analysis of organic connectivity.
        return Chem.RemoveHs(fragment)
    
    def carbon_count(fragment):
        # Count number of carbon atoms in a fragment (after removing H).
        frag = removeHs(fragment)
        return sum(1 for atom in frag.GetAtoms() if atom.GetSymbol() == "C")
    
    def fragment_is_organic(fragment):
        """
        Returns True if the fragment (after removing H) shows signs of extra organic features,
        e.g. contains any aromatic atom or a C=C double bond that is not simply part of a carbonyl.
        """
        frag = removeHs(fragment)
        # if no carbon, it's not organic.
        if not any(atom.GetSymbol() == "C" for atom in frag.GetAtoms()):
            return False
        # If any atom is aromatic, flag as organic.
        if any(atom.GetIsAromatic() for atom in frag.GetAtoms()):
            return True
        # Check each double bond between two carbons.
        for bond in frag.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                a1 = bond.GetBeginAtom()
                a2 = bond.GetEndAtom()
                if a1.GetSymbol() == "C" and a2.GetSymbol() == "C":
                    # If one of the carbons is in a carbonyl (C(=O)), then ignore.
                    def in_carbonyl(atom):
                        for nbr in atom.GetNeighbors():
                            b = frag.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                            if nbr.GetSymbol() == "O" and b and b.GetBondType() == Chem.BondType.DOUBLE:
                                return True
                        return False
                    if not (in_carbonyl(a1) or in_carbonyl(a2)):
                        return True
        return False

    def is_carboxylate(fragment):
        """
        Checks if a fragment contains a carboxylate substructure.
        Uses the SMARTS "[CX3](=O)[O-]".
        """
        carboxylate_smarts = Chem.MolFromSmarts("[CX3](=O)[O-]")
        return fragment.HasSubstructMatch(carboxylate_smarts)
    
    def allowed_organic_counterion(fragment, metal_is_alkali):
        """
        For a fragment that contains carbon but no allowed metal,
        decide if it is an allowed counterion.
         - If it is a carboxylate then for alkali metals only a one-carbon carboxylate (formate) is allowed;
           for other metals allow 1-3 carbons.
         - If it is not carboxylate, then if it is large (>=8 carbons) and shows no extra organic features, allow it.
        """
        cnt = carbon_count(fragment)
        # First require that it resembles a carboxylate substructure.
        if is_carboxylate(fragment):
            if metal_is_alkali:
                return cnt == 1
            else:
                return cnt in {1, 2, 3}
        # Not a carboxylate: if it is very large and not extra “organic”, we tolerate it.
        if cnt >= 8 and not fragment_is_organic(fragment):
            return True
        return False

    def is_metal_aqua_complex(fragment):
        """
        Returns True if the fragment is a metal-containing complex where the only non-metal atoms
        are hydrogen or oxygen. In our heuristic, these do not count as acceptable 
        (we want either an isolated metal ion or a salt with non-water counterions).
        """
        has_metal = False
        has_non_water_atoms = False
        for atom in fragment.GetAtoms():
            sym = atom.GetSymbol()
            if sym in allowed_metals:
                has_metal = True
            else:
                # If the atom is something other than H or O:
                if sym not in {"H", "O"}:
                    has_non_water_atoms = True
        # If we have a metal and no atoms besides H/O (and more than one atom overall) then it is water coordinated.
        if has_metal and (fragment.GetNumAtoms() > 1) and (not has_non_water_atoms):
            return True
        return False

    # ---------------- Main classification ----------------
    
    accepted_metal_msgs = []   # messages describing accepted metal ion fragments
    rejected_due_to_metal_complex = False
    # Process each fragment to find acceptable metal-containing fragments.
    for frag in frags:
        # Skip water fragments.
        if is_water(frag):
            continue

        frag_noHs = removeHs(frag)
        atoms = frag_noHs.GetAtoms()
        contains_allowed_metal = False
        for atom in atoms:
            sym = atom.GetSymbol()
            if sym not in allowed_metals:
                continue
            # If this fragment contains carbon it is too organic (unless it is an allowed counterion) 
            # so do not count metal ions embedded in an organic fragment.
            if carbon_count(frag) > 0:
                continue
            # Check for isolated metal ion: one atom with a positive formal charge.
            if frag.GetNumAtoms() == 1 and atom.GetFormalCharge() > 0:
                accepted_metal_msgs.append(f"Found isolated metal ion: {sym}{atom.GetFormalCharge()}")
                contains_allowed_metal = True
                break
            # Else if it is a multi‐atom inorganic fragment (no carbon) but with overall positive charge.
            overall_charge = Chem.GetFormalCharge(frag)
            if overall_charge > 0:
                accepted_metal_msgs.append(f"Found inorganic metal in non‐organic fragment: {sym}{'' if atom.GetFormalCharge()==0 else atom.GetFormalCharge()}")
                contains_allowed_metal = True
                break
        # Also check if the fragment is a metal–aqua complex. If so, disqualify.
        if contains_allowed_metal and is_metal_aqua_complex(frag):
            # Mark that we found a metal complex whose other atoms are only H/O (water ligands).
            rejected_due_to_metal_complex = True
        # Remember only fragments that successfully contained allowed metal and were NOT water complexes.
    if not accepted_metal_msgs:
        return False, "No acceptable metal nutrient ion found"
    if rejected_due_to_metal_complex:
        return False, "Metal nutrient ion found only as water‐coordinated (aqua) complex"

    # Determine if any accepted metal is alkali
    metal_is_alkali = any(msg.split()[3][0:-1] in alkali_metals for msg in accepted_metal_msgs)
    
    # Now process counterion (non–metal) fragments.
    counterion_errors = []
    for frag in frags:
        if is_water(frag):
            continue
        # If the fragment contains any allowed metal (and no carbon) we already processed it.
        # Otherwise, if the fragment has carbon then check if it is acceptable.
        frag_noHs = removeHs(frag)
        metals_in_frag = any(atom.GetSymbol() in allowed_metals for atom in frag_noHs.GetAtoms())
        if metals_in_frag:
            continue  # already handled as metal fragment
        if carbon_count(frag) > 0:
            # If the fragment appears overtly organic, disqualify it.
            if fragment_is_organic(frag):
                return False, "Disqualifying organic fragment detected: " + Chem.MolToSmiles(frag_noHs)
            if not allowed_organic_counterion(frag, metal_is_alkali):
                counterion_errors.append("Disallowed counterion fragment: " + Chem.MolToSmiles(frag_noHs))
    if counterion_errors:
        return False, ("Found metal nutrient ion(s): " + " ; ".join(accepted_metal_msgs)
                       + " but with counterion issue(s): " + " ; ".join(counterion_errors))
    
    return True, "Found metal nutrient ion(s): " + " ; ".join(accepted_metal_msgs)


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
      # Some false positives (expected to be rejected)
      "hexaaquamagnesium(2+)": "[H][O]([H])[Mg++]([O]([H])[H])([O]([H])[H])([O]([H])[H])([O]([H])[H])[O]([H])[H]",
      "Taurocholic acid sodium salt hydrate": "S([O-])(=O)(=O)CCNC(=O)CC[C@H]([C@@]1([C@@]2([C@]([C@]3([C@@]([C@@]4([C@](C[C@H]3O)(C[C@H](O)CC4)[H])C)(C[C@@H]2O)[H])[H])(CC1)[H])C)[H])C.[Na+].O",
      "sodium chlorite": "[Na+].[O-][Cl]=O",
      "sodium hydroxide": "[OH-].[Na+]",
      "sodium 3-aminopropyl 2-acetamido-2-deoxy-alpha-D-glucose-1-phosphate": "[Na+].CC(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]1OP([O-])(=O)OCCCN",
      "potassium tetrabromoaurate": "[K+].Br[Au-](Br)(Br)Br",
      "calcium titanate": "[Ca+2].[Ti+4].[O-2].[O-2].[O-2]",
      "aluminium phosphide": "[Al+3].[P-3]",
      # (More cases omitted for brevity.)
    }
    
    for name, smi in test_smiles.items():
        result, reason = is_mineral_nutrient(smi)
        print(f"{name}: {result}, {reason}")