"""
benchmarks/posada_2001_benchmark/scenarios.py
============================================
Defines the experimental grid replicating Posada & Crandall (2001, PNAS):
- Simulations I (Power Analysis): 4 theta levels x 5 rho levels = 20 conditions
- Simulations II (False Positive / Heterotachy): 4 theta levels x 3 alpha levels = 12 conditions
Total: 32 parameter configurations x 100 replicates = 3,200 simulation runs.
"""

from dataclasses import dataclass
from typing import Optional, List

@dataclass
class PosadaScenario:
    scenario_id: str
    category: str       # 'power' or 'false_positive'
    theta: float        # Population mutation parameter (10, 50, 100, 200)
    rho: float          # Population recombination parameter (0, 1, 4, 16, 64)
    alpha: Optional[float] # Gamma shape parameter (None for uniform, 2.0, 0.5, 0.05)
    n_taxa: int = 10
    length_nt: int = 1000
    n_reps: int = 100

def build_posada_scenarios() -> List[PosadaScenario]:
    scenarios = []
    
    # -------------------------------------------------------------------------
    # Simulations I: Power Grid (alpha = None, uniform rate across sites)
    # -------------------------------------------------------------------------
    for theta in [10.0, 50.0, 100.0, 200.0]:
        for rho in [0.0, 1.0, 4.0, 16.0, 64.0]:
            sid = f"power_theta{int(theta)}_rho{int(rho)}"
            scenarios.append(
                PosadaScenario(
                    scenario_id=sid,
                    category="power",
                    theta=theta,
                    rho=rho,
                    alpha=None,
                    n_taxa=10,
                    length_nt=1000,
                    n_reps=100
                )
            )
            
    # -------------------------------------------------------------------------
    # Simulations II: False Positive Grid (rho = 0.0, rate heterogeneity)
    # Note: alpha = None (uniform) is already included in power grid as rho=0.0
    # Here we add alpha in {2.0, 0.5, 0.05}
    # -------------------------------------------------------------------------
    for theta in [10.0, 50.0, 100.0, 200.0]:
        for alpha in [2.0, 0.5, 0.05]:
            alpha_str = f"{alpha}".replace(".", "p")
            sid = f"fp_theta{int(theta)}_alpha{alpha_str}"
            scenarios.append(
                PosadaScenario(
                    scenario_id=sid,
                    category="false_positive",
                    theta=theta,
                    rho=0.0,
                    alpha=alpha,
                    n_taxa=10,
                    length_nt=1000,
                    n_reps=100
                )
            )
            
    return scenarios
