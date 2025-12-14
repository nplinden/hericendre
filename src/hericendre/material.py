from pydantic import BaseModel

class Material(BaseModel):
    name: str = "Material"
    concentrations: dict[str, float] | None = None
    microxs: str | None = None

    def to_dict(self):
        data = {"name": self.name}
        if self.concentrations is not None:
            data["concentrations"] = self.concentrations
        if self.microxs is not None:
            data["microxs"] = self.microxs
        return data