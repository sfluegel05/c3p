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
       (bool, str): Tuple with overall classification and an explanation.
    """
    # Parse the SMILES string. If it fails, return an error.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # For fragment evaluation, work with fragments as separate molecules.
    # (Don’t add explicit hydrogens for the purpose of fragment splitting.)
    frags = Chem.rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    
    # Define allowed metal symbols (as strings) and which ones are considered alkali metals.
    allowed_metals = {"Pd", "K", "Fe", "Cs", "Ca", "Zn", "Al", "Mg", "Sb", "Ba", "Na", "La"}
    alkali_metals = {"Na", "K", "Cs"}
    
    # ---------------- Helper functions ----------------
    
    def is_water(fragment):
        """
        Returns True if the fragment corresponds to water.
        We use the molecular formula from rdMolDescriptors.
        """
        formula = rdMolDescriptors.CalcMolFormula(fragment)
        return formula == "H2O"
    
    def removeHs(fragment):
        """
        Returns a fragment with explicit hydrogens removed.
        """
        return Chem.RemoveHs(fragment)
    
    def carbon_count(fragment):
        """
        Count the number of carbon atoms in the fragment (after removing hydrogens).
        """
        frag = removeHs(fragment)
        return sum(1 for atom in frag.GetAtoms() if atom.GetSymbol() == "C")
    
    def fragment_is_organic(fragment):
        """
        Returns True if the fragment (after removing hydrogens) shows signs of having extra organic features,
        e.g. contains any aromatic atom or a C=C double bond that is not simply part of a carbonyl.
        """
        frag = removeHs(fragment)
        # If no carbon is present, we do not consider it organic.
        if not any(atom.GetSymbol() == "C" for atom in frag.GetAtoms()):
            return False
        # If any atom is aromatic, flag as organic.
        if any(atom.GetIsAromatic() for atom in frag.GetAtoms()):
            return True
        # Check for unconjugated carbon–carbon double bonds.
        for bond in frag.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                a1 = bond.GetBeginAtom()
                a2 = bond.GetEndAtom()
                if a1.GetSymbol() == "C" and a2.GetSymbol() == "C":
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
        Checks if a fragment contains a carboxylate substructure using the SMARTS "[CX3](=O)[O-]".
        """
        carboxylate_smarts = Chem.MolFromSmarts("[CX3](=O)[O-]")
        return fragment.HasSubstructMatch(carboxylate_smarts)
    
    def allowed_organic_counterion(fragment, metal_is_alkali):
        """
        For a fragment that contains carbon but no allowed metal, decide if it is an allowed counterion.
         - If it is carboxylate then for alkali metals only a one-carbon carboxylate (formate) is allowed;
           for other metals allow 1-3 carbons.
         - Otherwise, if the fragment is large (>=8 carbons) and shows no extra organic features, allow it.
        """
        cnt = carbon_count(fragment)
        if is_carboxylate(fragment):
            if metal_is_alkali:
                return cnt == 1
            else:
                return cnt in {1, 2, 3}
        if cnt >= 8 and not fragment_is_organic(fragment):
            return True
        return False

    def is_metal_aqua_complex(fragment):
        """
        Returns True if the fragment is a metal-containing complex where the non-metal atoms are only hydrogen or oxygen.
        Such water‐coordinated metal complexes are not acceptable.
        """
        has_metal = False
        has_non_water_atoms = False
        for atom in fragment.GetAtoms():
            sym = atom.GetSymbol()
            if sym in allowed_metals:
                has_metal = True
            else:
                if sym not in {"H", "O"}:
                    has_non_water_atoms = True
        # If the fragment has a metal and all its other atoms are H or O (and more than one atom overall), it is water‐coordinated.
        if has_metal and fragment.GetNumAtoms() > 1 and not has_non_water_atoms:
            return True
        return False

    # ---------------- Main classification ----------------
    
    accepted_metal_msgs = []    # messages describing acceptable metal fragments
    accepted_metal_symbols = [] # record the metal symbols found
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
            # If the fragment contains carbon, skip it (metal ions embedded in organic fragments are disqualified).
            if carbon_count(frag) > 0:
                continue
            # Check for isolated metal ion: one atom with a positive formal charge.
            if frag.GetNumAtoms() == 1 and atom.GetFormalCharge() > 0:
                accepted_metal_msgs.append(f"Found isolated metal ion: {sym}{atom.GetFormalCharge()}")
                accepted_metal_symbols.append(sym)
                contains_allowed_metal = True
                break
            # Otherwise, if it is a multi‐atom inorganic fragment (no carbon) with overall positive charge.
            overall_charge = Chem.GetFormalCharge(frag)
            if overall_charge > 0:
                accepted_metal_msgs.append(f"Found inorganic metal fragment: {sym}")
                accepted_metal_symbols.append(sym)
                contains_allowed_metal = True
                break
        # Also check if the fragment is a metal–aqua complex. If so, mark rejection.
        if contains_allowed_metal and is_metal_aqua_complex(frag):
            rejected_due_to_metal_complex = True
    if not accepted_metal_msgs:
        return False, "No acceptable metal nutrient ion found"
    if rejected_due_to_metal_complex:
        return False, "Metal nutrient ion found only as water‐coordinated (aqua) complex"
    
    # Determine if any accepted metal is alkali
    metal_is_alkali = any(m in alkali_metals for m in accepted_metal_symbols)
    
    # Process counterion (non‐metal) fragments.
    counterion_errors = []
    for frag in frags:
        if is_water(frag):
            continue
        frag_noHs = removeHs(frag)
        # Skip fragments that contain allowed metals (already processed) but check others.
        if any(atom.GetSymbol() in allowed_metals for atom in frag_noHs.GetAtoms()):
            continue
        if carbon_count(frag) > 0:
            # If the fragment appears overtly organic then disqualify.
            if fragment_is_organic(frag):
                return False, "Disqualifying organic fragment detected: " + Chem.MolToSmiles(frag_noHs)
            if not allowed_organic_counterion(frag, metal_is_alkali):
                counterion_errors.append("Disallowed counterion fragment: " + Chem.MolToSmiles(frag_noHs))
    if counterion_errors:
        error_msg = ("Found metal nutrient ion(s): " + " ; ".join(accepted_metal_msgs) +
                     " but with counterion issue(s): " + " ; ".join(counterion_errors))
        return False, error_msg
    
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
      "calcium hydrogenphosphate": "[Ca++].[H]OP([O-])([O-])=O"
    }
    
    for name, smi in test_smiles.items():
        result, reason = is_mineral_nutrient(smi)
        print(f"{name}: {result}, {reason}")