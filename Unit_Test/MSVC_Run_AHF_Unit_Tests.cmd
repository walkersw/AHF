rem script to run the unit tests.

@setlocal EnableDelayedExpansion

@for %%i in (test_vertex2halffacet_class.exe, test_2d_manifold_mesh_1.exe, test_2d_non_manifold_1.exe, test_2d_non_manifold_2.exe) do @(
  @echo --------------------------------------------------
  call %%i
  @if !errorlevel! equ 0 (
    @echo *Successfully* completed this unit test: %%i
  ) else (
    @echo The error code is:
	@echo !errorlevel!
    @echo This unit test *failed*: %%i
	@echo --------------------------------------------------
	exit /b
  )
)
@echo --------------------------------------------------

