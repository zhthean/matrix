import os

from conan import ConanFile


class MatrixRecipe(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeToolchain", "CMakeDeps"

    def requirements(self):
        self.requires("gtest/1.16.0")

    def layout(self):
        self.folders.generators = os.path.join("build", "generators")
        self.folders.build = os.path.join("build", str(self.settings.build_type))
