from dataclasses import dataclass


@dataclass(frozen=True)
class Material:
    """Isotropic benchmark material; values are inputs, not fitted data."""

    rho: float = 2700.0
    cL: float = 6320.0
    cT: float = 3130.0
    thickness: float = 1.0e-3

    @property
    def h(self) -> float:
        return self.thickness / 2

    @property
    def lame_mu(self) -> float:
        return self.rho * self.cT**2

    @property
    def lame_lambda(self) -> float:
        return self.rho * (self.cL**2 - 2 * self.cT**2)

    @property
    def poisson(self) -> float:
        return self.lame_lambda / (2 * (self.lame_lambda + self.lame_mu))

    @property
    def young(self) -> float:
        return self.lame_mu * (3 * self.lame_lambda + 2 * self.lame_mu) / (
            self.lame_lambda + self.lame_mu
        )

    @property
    def plate_speed(self) -> float:
        return (self.young / (self.rho * (1 - self.poisson**2))) ** 0.5
