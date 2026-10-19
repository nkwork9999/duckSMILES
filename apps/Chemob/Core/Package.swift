// swift-tools-version: 6.0
import PackageDescription
let package = Package(name: "ChemobCore", platforms: [.iOS(.v17), .macOS(.v14)], products: [.library(name: "ChemobCore", targets: ["ChemobCore"])], targets: [.target(name: "ChemobCore"), .testTarget(name: "ChemobCoreTests", dependencies: ["ChemobCore"])])
