import random
from typing import List, Optional
from biobridge.definitions.cells.muscle_cell import MuscleCell
from biobridge.blocks.tissue import Tissue


class MuscleTissue(Tissue):
    def __init__(self, name: str, cells: Optional[List[MuscleCell]] = None, cancer_risk: float = 0.001):
        """
        Initialize a new MuscleTissue object, inheriting from Tissue.

        :param name: Name of the muscle tissue
        :param cells: List of MuscleCell objects that make up the muscle tissue
        :param cancer_risk: The risk of cancer for this tissue type
        """
        super().__init__(name, tissue_type="muscle", cells=cells, cancer_risk=cancer_risk)
        self.contraction_rate = 0.3  # Percentage of muscle cells that contract simultaneously

    def contract_muscle(self):
        """Simulate the contraction of the muscle tissue."""
        num_cells = self.get_cell_count()
        contracting_cells = int(num_cells * self.contraction_rate)
        
        # Randomly select muscle cells to contract
        for cell in random.sample(self.cells, contracting_cells):
            if isinstance(cell, MuscleCell) and not cell.contracted:
                cell.contract(1)

    def relax_muscle(self):
        """Simulate the relaxation of the muscle tissue."""
        for cell in self.cells:
            if isinstance(cell, MuscleCell) and cell.contracted:
                cell.relax()

    def simulate_contraction_cycle(self):
        """Simulate a full cycle of muscle contraction and relaxation."""
        print(f"Simulating contraction for {self.name}...")
        self.contract_muscle()
        self.tissue_metabolism()  # During contraction, the muscle cells metabolize
        self.relax_muscle()
        
    def simulate_time_step(self, external_factors: List[tuple] = None):
        """Override the simulate_time_step method to include muscle contraction."""
        super().simulate_time_step(external_factors)
        # Add muscle-specific behavior
        if random.random() < 0.5:  # 50% chance of simulating a contraction
            self.simulate_contraction_cycle()

    def get_total_contractile_force(self) -> float:
        """Calculate and return the total contractile force generated by the muscle tissue."""
        return sum(cell.contractile_force for cell in self.cells if isinstance(cell, MuscleCell))

    def describe(self) -> str:
        """Provide a detailed description of the muscle tissue."""
        base_description = super().describe()
        total_force = self.get_total_contractile_force()
        return f"{base_description}\nTotal Contractile Force: {total_force:.2f}"