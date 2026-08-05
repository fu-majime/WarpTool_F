[CmdletBinding()]
param(
    [ValidateSet("Debug", "Release")]
    [string]$Configuration = "Release"
)

$ErrorActionPreference = "Stop"
$repositoryRoot = $PSScriptRoot
$profile = $Configuration.ToLowerInvariant()
$cargoArguments = @("build", "--package", "puppet_mesh")
if ($Configuration -eq "Release") {
    $cargoArguments += "--release"
}

Push-Location -LiteralPath $repositoryRoot
try {
    & cargo @cargoArguments
    if ($LASTEXITCODE -ne 0) {
        throw "cargo build failed with exit code $LASTEXITCODE"
    }

    $distributionDirectory = Join-Path $repositoryRoot "dist"
    New-Item -ItemType Directory -Path $distributionDirectory -Force | Out-Null

    $moduleSource = Join-Path $repositoryRoot "target\$profile\puppet_mesh.dll"
    $moduleDestination = Join-Path $distributionDirectory "puppet_mesh.mod2"
    $scriptSource = Join-Path $repositoryRoot "lua\@WarpTool_F.anm2"
    $waveWarpSource = Join-Path $repositoryRoot "lua\WaveWarp_F.lua"

    Copy-Item -LiteralPath $moduleSource -Destination $moduleDestination -Force
    Copy-Item -LiteralPath $scriptSource -Destination $distributionDirectory -Force
    Copy-Item -LiteralPath $waveWarpSource -Destination $distributionDirectory -Force

    Write-Host "Created AviUtl2 package: $distributionDirectory"
}
finally {
    Pop-Location
}
