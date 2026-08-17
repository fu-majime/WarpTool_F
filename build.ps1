[CmdletBinding()]
param(
    [ValidateSet("Debug", "Release")]
    [string]$Configuration = "Release"
)

$ErrorActionPreference = "Stop"
$repositoryRoot = $PSScriptRoot
$profile = $Configuration.ToLowerInvariant()
$cargoArguments = @(
    "build",
    "--package", "puppet_geometry",
    "--package", "puppet_tool_aux",
    "--package", "script_edit_bridge"
)
if ($Configuration -eq "Release") {
    $cargoArguments += "--release"
}

Push-Location -LiteralPath $repositoryRoot
try {
    & aulua build
    if ($LASTEXITCODE -ne 0) {
        throw "aulua build failed with exit code $LASTEXITCODE"
    }

    & cargo @cargoArguments
    if ($LASTEXITCODE -ne 0) {
        throw "cargo build failed with exit code $LASTEXITCODE"
    }

    $tempDistDirectory = Join-Path $repositoryRoot "target\dist"
    $pluginDirectory = Join-Path $tempDistDirectory "Plugin\WarpTool_F"
    $scriptDirectory = Join-Path $tempDistDirectory "Script\WarpTool_F"

    if (Test-Path -LiteralPath $tempDistDirectory) {
        Remove-Item -LiteralPath $tempDistDirectory -Recurse -Force | Out-Null
    }
    New-Item -ItemType Directory -Path $pluginDirectory -Force | Out-Null
    New-Item -ItemType Directory -Path $scriptDirectory -Force | Out-Null

    $moduleSource = Join-Path $repositoryRoot "target\$profile\puppet_geometry.dll"
    $moduleDestination = Join-Path $scriptDirectory "puppet_geometry.mod2"
    $bridgeSource = Join-Path $repositoryRoot "target\$profile\script_edit_bridge.dll"
    $bridgeDestination = Join-Path $pluginDirectory "ScriptEditBridge.aux2"
    $puppetEditorSource = Join-Path $repositoryRoot "target\$profile\puppet_tool_aux.dll"
    $puppetEditorDestination = Join-Path $pluginDirectory "WarpTool_F.aux2"
    $scriptSource = Join-Path $repositoryRoot "target\lua\@WarpTool_F.anm2"
    $waveWarpSource = Join-Path $repositoryRoot "target\lua\WaveWarpGeometry.lua"

    Copy-Item -LiteralPath $moduleSource -Destination $moduleDestination -Force
    Copy-Item -LiteralPath $bridgeSource -Destination $bridgeDestination -Force
    Copy-Item -LiteralPath $puppetEditorSource -Destination $puppetEditorDestination -Force
    Copy-Item -LiteralPath $scriptSource -Destination $scriptDirectory -Force
    Copy-Item -LiteralPath $waveWarpSource -Destination $scriptDirectory -Force

    $distributionDirectory = Join-Path $repositoryRoot "dist"
    if (Test-Path -LiteralPath $distributionDirectory) {
        Remove-Item -LiteralPath $distributionDirectory -Recurse -Force | Out-Null
    }
    New-Item -ItemType Directory -Path $distributionDirectory -Force | Out-Null

    $zipPath = Join-Path $distributionDirectory "WarpTool_F.au2pkg.zip"
    Add-Type -AssemblyName System.IO.Compression.FileSystem
    [System.IO.Compression.ZipFile]::CreateFromDirectory($tempDistDirectory, $zipPath)

    Remove-Item -LiteralPath $tempDistDirectory -Recurse -Force | Out-Null

    Write-Host "Created AviUtl2 package zip: $zipPath"
}
finally {
    Pop-Location
}
