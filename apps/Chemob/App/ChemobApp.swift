import SwiftUI
import UIKit
import UniformTypeIdentifiers
import ChemobCore

enum ChemobStyle {
    static let ink = Color(red: 0.07, green: 0.12, blue: 0.13)
    static let paper = Color(red: 0.96, green: 0.97, blue: 0.94)
    static let lime = Color(red: 0.78, green: 0.94, blue: 0.28)
    static let muted = Color(red: 0.35, green: 0.40, blue: 0.39)
    static let warning = Color(red: 0.62, green: 0.25, blue: 0.03)
    static let oxygen = Color(red: 0.72, green: 0.17, blue: 0.12)
}

@main struct ChemobApp: App {
    var body: some Scene { WindowGroup { WorkbenchView() } }
}

struct WorkbenchView: View {
    @StateObject private var model = Workbench()
    @State private var paste = false
    @State private var importing = false
    @State private var importChoice = false
    @State private var fileReplacing = false
    @State private var filtersVisible = false
    @State private var exportURL: URL?
    @State private var selected: CompoundInput?
    @Environment(\.dynamicTypeSize) private var textSize
    @ScaledMetric(relativeTo: .largeTitle) private var heroTitleSize = 39
    var body: some View {
        NavigationStack {
            Group {
                if textSize.isAccessibilitySize { accessibleContent } else {
            GeometryReader { geometry in
                let wide = geometry.size.width >= 700
                if wide && !model.inputs.isEmpty {
                    InlineWorkbench(model: model, add: { paste = true }, open: { requestFileImport() }, export: beginExport).padding(20)
                } else {
                VStack(spacing: 0) {
                    VStack(alignment: .leading, spacing: 12) {
                        HStack {
                            HStack(spacing: 8) {
                                Image(systemName: "hexagon.fill").foregroundStyle(ChemobStyle.ink)
                                Text("Chemob").font(.system(.title2, design: .rounded).weight(.heavy))
                            }
                            Spacer()
                            Label("On device", systemImage: "checkmark.shield")
                                .font(.caption.weight(.medium)).foregroundStyle(ChemobStyle.muted)
                        }
                        if !model.inputs.isEmpty {
                            HStack(alignment: .firstTextBaseline) {
                                Text("Compounds").font(.system(.title, design: .rounded).bold())
                                Text(model.inputs.count, format: .number).font(.system(.title, design: .monospaced)).foregroundStyle(ChemobStyle.muted)
                                Spacer()
                                if model.busy { ProgressView().tint(ChemobStyle.ink) }
                                else if let elapsed = model.elapsed {
                                    Text(elapsed, format: .number.precision(.fractionLength(2))) + Text(" s")
                                }
                            }.monospacedDigit()
                        }
                        HStack(spacing: 12) {
                            Button { paste = true } label: {
                                Label("Add SMILES", systemImage: "plus").font(.subheadline.bold()).frame(maxWidth: .infinity, minHeight: 48)
                            }.buttonStyle(WorkbenchButtonStyle(primary: true)).accessibilityIdentifier("pasteSMILES")
                            Button { requestFileImport() } label: {
                                Label("Import", systemImage: "folder").font(.subheadline.bold()).frame(minWidth: 96, minHeight: 48)
                            }.buttonStyle(WorkbenchButtonStyle(primary: false)).accessibilityIdentifier("importFile")
                        }
                        if !model.inputs.isEmpty {
                            HStack(spacing: 10) {
                                Image(systemName: "magnifyingglass").accessibilityHidden(true)
                                TextField("Search name or SMILES", text: $model.query)
                                    .textInputAutocapitalization(.never).autocorrectionDisabled().accessibilityIdentifier("compoundSearch")
                                if !model.query.isEmpty {
                                    Button { model.query = "" } label: { Image(systemName: "xmark.circle.fill").frame(width: 44, height: 44) }
                                        .accessibilityLabel("Clear search")
                                }
                                NumericFilterToggle(active: model.numericFilter != nil, expanded: $filtersVisible)
                            }.padding(.horizontal, 14).frame(minHeight: 48).background(.white, in: RoundedRectangle(cornerRadius: 12))
                            if filtersVisible { NumericFilterEditor(model: model, compact: !wide) { filtersVisible = false } }
                            else if model.numericFilter != nil { NumericFilterSummary(model: model) }
                        }
                    }.padding(.horizontal, wide ? 28 : 20).padding(.vertical, 16)
                    if model.inputs.isEmpty {
                        ScrollView {
                            VStack(alignment: .leading, spacing: 24) {
                                hero
                                Button { model.examples() } label: {
                                    HStack {
                                        Image(systemName: "arrow.up.right.square")
                                        VStack(alignment: .leading, spacing: 4) {
                                            Text("Try a sample set").font(.subheadline.bold())
                                            Text("Six inputs, including one invalid SMILES").font(.caption).foregroundStyle(ChemobStyle.muted)
                                        }
                                        Spacer()
                                        Image(systemName: "arrow.right")
                                    }.padding(18)
                                }.buttonStyle(WorkbenchButtonStyle(primary: false)).accessibilityIdentifier("loadExamples")
                            }.padding(24).frame(maxWidth: 640).frame(maxWidth: .infinity)
                        }
                    } else {
                        HStack {
                            VStack(alignment: .leading, spacing: 3) {
                                if model.busy { Text("Processed \(model.results.count) of \(model.inputs.count)") }
                                else if model.exporting { Text("Exporting CSV") }
                                else if model.reusedCount > 0 { Text("Reused \(model.reusedCount) results") }
                                else { Text("Calculated locally") }
                                if !model.query.isEmpty || model.numericFilter != nil {
                                    Group { if model.visible.count == 1 { Text("1 match") } else { Text("\(model.visible.count) matches") } }.foregroundStyle(ChemobStyle.muted)
                                }
                            }.font(.caption).accessibilityIdentifier("calculationStatus")
                            Spacer()
                            if model.busy { Button("Cancel") { model.cancel() }.frame(minHeight: 44) }
                            else {
                                CSVExportControl(model: model, action: beginExport)
                            }
                        }.padding(.horizontal, wide ? 28 : 20)
                        HStack(spacing: 16) {
                            Text("COMPOUND").frame(maxWidth: .infinity, alignment: .leading)
                            Text("MW").frame(width: 72, alignment: .trailing)
                            if wide {
                                Text("cLogP").frame(width: 76, alignment: .trailing)
                                Text("TPSA / Å²").frame(width: 90, alignment: .trailing)
                                Text("HBD / HBA").frame(width: 90, alignment: .trailing)
                            }
                            Color.clear.frame(width: 12, height: 1)
                        }.font(.caption2.bold()).foregroundStyle(ChemobStyle.muted).padding(.horizontal, wide ? 44 : 36).padding(.vertical, 12)
                        ScrollView {
                            LazyVStack(spacing: 0) {
                                ForEach(model.visible) { input in
                                    Button { selected = input } label: { CompoundCard(input: input, result: model.result(for: input), wide: wide) }
                                        .buttonStyle(CompoundRowStyle()).accessibilityIdentifier("compound-\(input.id)")
                                }
                            }.background(.white, in: RoundedRectangle(cornerRadius: 16))
                                .padding(.horizontal, wide ? 28 : 20).padding(.bottom, 20)
                            if model.visible.isEmpty {
                                ContentUnavailableView.search(text: model.query)
                            }
                        }
                    }
                }.frame(maxWidth: 1200).frame(maxWidth: .infinity)
                }
            }
                }
            }
            .background(ChemobStyle.paper).foregroundStyle(ChemobStyle.ink)
            .toolbar(.hidden, for: .navigationBar)
            .sheet(isPresented: $paste) { PasteView(existingCount: model.inputs.count) { text, format, replacing in model.importText(text, format: format, replacing: replacing) } }
            .confirmationDialog("Import compounds", isPresented: $importChoice, titleVisibility: .visible) {
                Button("Add to library") { fileReplacing = false; importing = true }
                Button("Replace library") { fileReplacing = true; importing = true }
                Button("Cancel", role: .cancel) { }
            }
            .sheet(item: $selected) { input in CompoundDetail(input: input, result: model.result(for: input)) }
            .sheet(isPresented: Binding(get: { exportURL != nil }, set: { if !$0 {
                if let exportURL { model.discardCSV(exportURL) }
                exportURL = nil
            } })) { if let exportURL { ActivitySheet(url: exportURL) } }
            .fileImporter(isPresented: $importing, allowedContentTypes: [.commaSeparatedText, .tabSeparatedText, .plainText, .data]) { response in
                do {
                    let url = try response.get(); let granted = url.startAccessingSecurityScopedResource()
                    defer { if granted { url.stopAccessingSecurityScopedResource() } }
                    let data = try Data(contentsOf: url, options: .mappedIfSafe)
                    guard data.count <= CompoundImport.maximumInputBytes, let text = String(data: data, encoding: .utf8) else { throw ImportIssue.tooLarge }
                    model.importText(text, format: url.pathExtension, replacing: fileReplacing)
                } catch { model.error = error.localizedDescription }
            }
            .alert("Check your input", isPresented: Binding(get: { model.error != nil }, set: { if !$0 { model.error = nil } })) { Button("OK") { model.error = nil } } message: { Text(model.error ?? "") }
            .task { model.restore() }
        }.tint(ChemobStyle.ink).preferredColorScheme(.light)
    }
    private func requestFileImport() {
        if model.inputs.isEmpty { fileReplacing = false; importing = true }
        else { importChoice = true }
    }
    private func beginExport() {
        Task {
            do { exportURL = try await model.exportCSV() }
            catch is CancellationError { }
            catch { model.error = error.localizedDescription }
        }
    }
    private var accessibleContent: some View {
        ScrollView {
            LazyVStack(alignment: .leading, spacing: 20) {
                Text("Chemob").font(.title2.bold())
                Label("On device", systemImage: "checkmark.shield").font(.caption)
                if model.inputs.isEmpty { hero } else {
                    Text("Compounds").font(.title2.bold())
                    Text(model.inputs.count, format: .number).font(.title3.monospacedDigit())
                }
                Button { paste = true } label: {
                    Label("Add SMILES", systemImage: "plus").frame(maxWidth: .infinity, alignment: .leading).padding(12)
                }.buttonStyle(WorkbenchButtonStyle(primary: true)).accessibilityIdentifier("pasteSMILES")
                Button { requestFileImport() } label: {
                    Label("Import", systemImage: "folder").frame(maxWidth: .infinity, alignment: .leading).padding(12)
                }.buttonStyle(WorkbenchButtonStyle(primary: false)).accessibilityIdentifier("importFile")
                if model.inputs.isEmpty {
                    Button("Try a sample set") { model.examples() }.padding(.vertical, 12).accessibilityIdentifier("loadExamples")
                } else {
                    TextField("Search name or SMILES", text: $model.query).textInputAutocapitalization(.never).autocorrectionDisabled()
                        .padding(12).background(.white, in: RoundedRectangle(cornerRadius: 12)).accessibilityIdentifier("compoundSearch")
                    if !model.query.isEmpty { Button("Clear search") { model.query = "" }.frame(minHeight: 44) }
                    NumericFilterToggle(active: model.numericFilter != nil, expanded: $filtersVisible)
                    if filtersVisible { NumericFilterEditor(model: model, compact: true) { filtersVisible = false } }
                    else if model.numericFilter != nil { NumericFilterSummary(model: model) }
                    if model.busy {
                        Text("Processed \(model.results.count) of \(model.inputs.count)")
                        Button("Cancel") { model.cancel() }.frame(minHeight: 44)
                    } else {
                        if model.exporting { Text("Exporting CSV") }
                        else if model.reusedCount > 0 { Text("Reused \(model.reusedCount) results") }
                        CSVExportControl(model: model, action: beginExport)
                    }
                    ForEach(model.visible) { input in
                        Button { selected = input } label: { CompoundCard(input: input, result: model.result(for: input)) }
                            .buttonStyle(CompoundRowStyle()).background(.white, in: RoundedRectangle(cornerRadius: 16)).accessibilityIdentifier("compound-\(input.id)")
                    }
                }
            }.padding(20).frame(maxWidth: 800).frame(maxWidth: .infinity)
        }
    }
    private var hero: some View {
        VStack(alignment: .leading, spacing: 12) {
            SpriteIntro().frame(height: 190).frame(maxWidth: .infinity)
            Text("Small molecules.\nBig momentum.").font(.system(size: heroTitleSize, weight: .bold, design: .rounded)).tracking(-1.5).fixedSize(horizontal: false, vertical: true)
            Text("Paste a batch. Compare properties.\nKeep your chemistry moving, offline.").font(.subheadline).foregroundStyle(ChemobStyle.muted).lineSpacing(5)
            LazyVGrid(columns: [GridItem(.adaptive(minimum: 140), alignment: .leading)], alignment: .leading, spacing: 8) { ForEach([String.localizedStringWithFormat(NSLocalizedString("%lld descriptors", comment: "Descriptor count"), Chemistry.descriptorCatalog.count), "Batch analysis", "CSV export"], id: \.self) { Text(LocalizedStringKey($0)).font(.caption.bold()).padding(.horizontal, 12).padding(.vertical, 8).background(ChemobStyle.lime, in: Capsule()) } }.padding(.top, 6)
        }
    }
}

struct WorkbenchButtonStyle: ButtonStyle {
    var primary: Bool
    func makeBody(configuration: Configuration) -> some View {
        configuration.label.foregroundStyle(ChemobStyle.ink)
            .background(primary ? ChemobStyle.lime : Color.white, in: RoundedRectangle(cornerRadius: 12))
            .opacity(configuration.isPressed ? 0.65 : 1)
    }
}

struct CompoundRowStyle: ButtonStyle {
    func makeBody(configuration: Configuration) -> some View {
        configuration.label.background(configuration.isPressed ? ChemobStyle.lime.opacity(0.25) : Color.clear)
    }
}

struct CompoundCard: View {
    let input: CompoundInput
    let result: AnalysisRow?
    var wide = false
    @Environment(\.dynamicTypeSize) private var textSize
    @Environment(\.locale) private var locale
    var body: some View {
        VStack(alignment: .leading, spacing: 8) {
            HStack(alignment: .center, spacing: 16) {
                VStack(alignment: .leading, spacing: 5) {
                    Text(input.name).font(.subheadline.weight(.semibold)).lineLimit(textSize.isAccessibilitySize ? nil : 1)
                    Text(result?.formula ?? input.smiles).font(.system(.caption, design: .monospaced)).foregroundStyle(ChemobStyle.muted).lineLimit(1)
                }.frame(maxWidth: .infinity, alignment: .leading)
                if !textSize.isAccessibilitySize {
                    value(result?.mw, width: 72)
                    if wide {
                        value(result?.logp, width: 76)
                        value(result?.tpsa, width: 90)
                        Text(result?.hbd.map { "\($0) / \(result?.hba ?? 0)" } ?? "—").font(.system(.subheadline, design: .monospaced)).frame(width: 90, alignment: .trailing)
                    }
                }
                Image(systemName: "chevron.right").font(.caption2.bold()).foregroundStyle(ChemobStyle.muted).frame(width: 12).accessibilityHidden(true)
            }
            if textSize.isAccessibilitySize {
                Text("MW: \(DescriptorFormatting.text(result?.mw, digits: 2, locale: locale)) g/mol").font(.body.monospacedDigit())
            }
            if result?.isValid == false {
                Label("Invalid or unsupported SMILES", systemImage: "exclamationmark.triangle").font(.caption).foregroundStyle(ChemobStyle.warning)
            } else if !wide, let result, result.isValid {
                (textSize.isAccessibilitySize ? AnyLayout(VStackLayout(alignment: .leading, spacing: 8)) : AnyLayout(HStackLayout(spacing: 16))) {
                    Text("cLogP \(DescriptorFormatting.text(result.logp, digits: 2, locale: locale))")
                    Text("TPSA \(DescriptorFormatting.text(result.tpsa, digits: 1, locale: locale)) Å²")
                }.font(.system(.caption2, design: .monospaced)).foregroundStyle(ChemobStyle.muted)
            }
        }.padding(.horizontal, 16).padding(.vertical, 12).contentShape(Rectangle())
            .overlay(alignment: .bottom) { Rectangle().fill(ChemobStyle.paper).frame(height: 1).padding(.horizontal, 16) }
    }
    private func value(_ number: Double?, width: CGFloat) -> some View {
        Text(DescriptorFormatting.text(number, digits: 2, locale: locale))
            .font(.system(.subheadline, design: .monospaced).weight(.medium)).frame(width: width, alignment: .trailing)
    }
}

struct PasteView: View {
    @Environment(\.dismiss) private var dismiss
    @State private var text = ""
    @State private var format = "smi"
    @State private var replacing = false
    let existingCount: Int
    let submit: (String, String, Bool) -> Void
    var body: some View {
        NavigationStack {
            VStack(alignment: .leading, spacing: 18) {
                Picker("Format", selection: $format) {
                    Text("SMILES").tag("smi")
                    Text("CSV").tag("csv")
                    Text("TSV").tag("tsv")
                }.pickerStyle(.segmented).accessibilityIdentifier("inputFormat")
                if existingCount > 0 {
                    Picker("Library", selection: $replacing) {
                        Text("Add to library").tag(false)
                        Text("Replace library").tag(true)
                    }.pickerStyle(.segmented).accessibilityIdentifier("inputMode")
                    if replacing { Text("The current list will be replaced.").font(.caption).foregroundStyle(ChemobStyle.warning) }
                }
                Text(LocalizedStringKey(format == "smi" ? "One SMILES per line, followed by an optional name." : "Use a smiles column and an optional name column.")).foregroundStyle(.secondary)
                Text(format == "smi" ? "CCO Ethanol\nc1ccccc1 Benzene" : format == "csv" ? "smiles,name\nCCO,Ethanol" : "smiles\tname\nCCO\tEthanol").font(.system(.callout, design: .monospaced)).foregroundStyle(.secondary)
                TextEditor(text: $text).font(.system(.body, design: .monospaced)).keyboardType(.asciiCapable).textInputAutocapitalization(.never).autocorrectionDisabled().padding(8).background(ChemobStyle.paper, in: RoundedRectangle(cornerRadius: 16)).accessibilityIdentifier("smilesEditor")
                Button { submit(text, format, replacing); dismiss() } label: { Text(LocalizedStringKey(replacing ? "Replace and analyze" : "Analyze batch")).font(.headline).frame(maxWidth: .infinity).padding(18) }.buttonStyle(.plain).foregroundStyle(.white).background(ChemobStyle.ink, in: RoundedRectangle(cornerRadius: 16)).disabled(text.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty).accessibilityIdentifier("analyzeInput")
            }.padding(24).navigationTitle(Text(LocalizedStringKey(format == "smi" ? "Add SMILES" : "Add compounds"))).navigationBarTitleDisplayMode(.inline).toolbar { ToolbarItem(placement: .cancellationAction) { Button("Close") { dismiss() } } }
        }
    }
}

struct CompoundDetail: View {
    @Environment(\.dismiss) private var dismiss
    @Environment(\.locale) private var locale
    let input: CompoundInput
    let result: AnalysisRow?
    @State private var graph: ChemicalDepiction?
    var body: some View {
        NavigationStack {
            ScrollView {
                VStack(alignment: .leading, spacing: 24) {
                    VStack(spacing: 10) {
                        if let graph { ChemicalStructurePreview(drawing: graph).frame(height: 220) }
                        Text(result?.formula ?? (result == nil ? String(localized: "Not calculated") : String(localized: "Check the structure"))).font(.system(.largeTitle, design: .rounded).bold())
                        Text(LocalizedStringKey(graph?.stereochemistryOmitted == true ? "2D structure · some stereochemistry not shown" : graph == nil ? "Preview unavailable" : "2D structure")).font(.caption).foregroundStyle(ChemobStyle.muted)
                    }.frame(maxWidth: .infinity).padding(24).background(.white, in: RoundedRectangle(cornerRadius: 26))
                    if let result, result.isValid {
                        DescriptorBrowser(row: result, columns: 2)
                    } else if result == nil {
                        Label("This molecule has not been calculated yet.", systemImage: "hourglass").foregroundStyle(ChemobStyle.muted)
                    } else { Label("Invalid or unsupported structure. The original SMILES is preserved.", systemImage: "exclamationmark.triangle").foregroundStyle(ChemobStyle.warning) }
                    Text("SMILES").font(.caption.bold()).foregroundStyle(ChemobStyle.muted)
                    Text(input.smiles).font(.system(.callout, design: .monospaced)).textSelection(.enabled)
                }.padding(24).frame(maxWidth: 720).frame(maxWidth: .infinity)
            }.background(ChemobStyle.paper).navigationTitle(input.name).navigationBarTitleDisplayMode(.inline).toolbar { ToolbarItem(placement: .confirmationAction) { Button("Done") { dismiss() } } }
                .task { let value = await ChemistryStore.shared.drawing(input.smiles); if !Task.isCancelled { graph = value } }
        }.tint(ChemobStyle.ink)
    }
    private func metric(_ name: String, _ value: Double?, _ unit: String, digits: Int = 2) -> some View {
        VStack(alignment: .leading, spacing: 8) {
            Text(LocalizedStringKey(name)).font(.caption).foregroundStyle(ChemobStyle.muted)
            Text(DescriptorFormatting.text(value, digits: digits, locale: locale)).font(.system(.title, design: .rounded).weight(.semibold)).minimumScaleFactor(0.65).lineLimit(1)
            Text(LocalizedStringKey(unit)).font(.caption2).foregroundStyle(ChemobStyle.muted)
        }.frame(maxWidth: .infinity, alignment: .leading).padding(18).background(.white, in: RoundedRectangle(cornerRadius: 20))
    }
}

struct ChemicalStructurePreview: View {
    let drawing: ChemicalDepiction
    var body: some View {
        Canvas { context, size in
            let minX=drawing.points.map(\.x).min() ?? 0, maxX=drawing.points.map(\.x).max() ?? 0
            let minY=drawing.points.map(\.y).min() ?? 0, maxY=drawing.points.map(\.y).max() ?? 0
            let scale=min(80, min(size.width / (maxX-minX+1.5), size.height / (maxY-minY+1.5)))
            func point(_ p: DrawingPoint) -> CGPoint {
                let screen=drawing.screenPoint(p)
                return CGPoint(x: size.width/2+(screen.x-(minX+maxX)/2)*scale,y: size.height/2+(screen.y+(minY+maxY)/2)*scale)
            }
            for stroke in drawing.strokes {
                var path=Path(); path.move(to: point(stroke.a)); path.addLine(to: point(stroke.b))
                context.stroke(path,with: .color(ChemobStyle.ink),style: StrokeStyle(lineWidth: max(1,scale*0.035),lineCap: .round,dash: stroke.dashed ? [scale*0.09,scale*0.06] : []))
            }
            for wedge in drawing.wedges {
                let d=wedge.b-wedge.a,normal=DrawingPoint(x: -d.y,y: d.x)*(0.15/d.length)
                if wedge.hashed {
                    for i in 1...7 {
                        let t=Double(i)/8,p=wedge.a+d*t
                        var path=Path();path.move(to: point(p-normal*t));path.addLine(to: point(p+normal*t))
                        context.stroke(path,with: .color(ChemobStyle.ink),lineWidth: max(1,scale*0.025))
                    }
                } else {
                    var path=Path();path.move(to: point(wedge.a));path.addLine(to: point(wedge.b-normal));path.addLine(to: point(wedge.b+normal));path.closeSubpath()
                    context.fill(path,with: .color(ChemobStyle.ink))
                }
            }
            let fontSize=max(10,min(19,scale*0.29))
            // Resolve actual glyph dimensions before masking bonds behind labels.
            for label in drawing.labels {
                let text=Text(label.text).font(.system(size: fontSize,weight: .medium))
                let resolved=context.resolve(text)
                let bounds=resolved.measure(in: CGSize(width: 500,height: 100))
                let center=point(label.point)
                let rect=CGRect(x: center.x-bounds.width/2-2,y: center.y-bounds.height/2-1,width: bounds.width+4,height: bounds.height+2)
                context.fill(Path(roundedRect: rect,cornerRadius: 2),with: .color(.white))
            }
            for label in drawing.labels {
                context.draw(Text(label.text).font(.system(size: fontSize,weight: .medium)).foregroundStyle(label.element == "O" ? ChemobStyle.oxygen : ChemobStyle.ink),at: point(label.point))
            }
        }.accessibilityLabel(Text(drawing.stereochemistryOmitted ? "2D structure · some stereochemistry not shown" : "2D structure"))
            .accessibilityIdentifier("chemicalDepiction")
    }
}

private struct CSVExportControl: View {
    @ObservedObject var model: Workbench
    let action: () -> Void
    var iconOnly = false
    var body: some View {
        Group {
            if model.exporting {
                Button { model.cancelExport() } label: {
                    HStack(spacing: 6) {
                        ProgressView().controlSize(.small)
                        if iconOnly { Image(systemName: "xmark") }
                        else { Text("Cancel") }
                    }.frame(width: iconOnly ? 44 : nil, height: 44)
                }.accessibilityLabel("Cancel CSV export")
                    .accessibilityValue(Text("Exporting CSV")).accessibilityIdentifier("cancelExport")
            } else {
                Button(action: action) {
                    if iconOnly { Image(systemName: "square.and.arrow.up").frame(width: 44, height: 44) }
                    else { Label("Export CSV", systemImage: "square.and.arrow.up").frame(minHeight: 44) }
                }.accessibilityLabel("Export CSV").accessibilityIdentifier("exportCSV")
            }
        }.font(.caption.bold()).buttonStyle(.plain)
    }
}

struct ActivitySheet: UIViewControllerRepresentable {
    let url: URL
    func makeUIViewController(context: Context) -> UIActivityViewController {
        let controller = UIActivityViewController(activityItems: [url], applicationActivities: nil)
        controller.view.accessibilityIdentifier = "csvShareSheet"
        return controller
    }
    func updateUIViewController(_ controller: UIActivityViewController, context: Context) { }
}

private struct NumericFilterToggle: View {
    let active: Bool
    @Binding var expanded: Bool
    var body: some View {
        Button { expanded.toggle() } label: {
            Image(systemName: active ? "line.3.horizontal.decrease.circle.fill" : "line.3.horizontal.decrease.circle")
                .font(.title3).frame(width: 44, height: 44)
                .background(active ? ChemobStyle.lime : .clear, in: RoundedRectangle(cornerRadius: 10))
        }.buttonStyle(.plain).accessibilityLabel("Numeric filter")
            .accessibilityValue(Text(LocalizedStringKey(active ? "Applied" : "Not applied")))
            .accessibilityIdentifier("numericFilter")
    }
}

private struct NumericFilterSummary: View {
    @ObservedObject var model: Workbench
    @Environment(\.locale) private var locale
    var body: some View {
        if let filter = model.numericFilter {
            let spec = Chemistry.descriptorCatalog.first { $0.id == filter.descriptorID }
            HStack(spacing: 6) {
                Text(LocalizedStringKey(spec?.title ?? filter.descriptorID)).lineLimit(1)
                Text(range(filter)).font(.caption.monospacedDigit()).lineLimit(1)
                if let unit = spec?.unit, !unit.isEmpty { Text(unit).lineLimit(1) }
                Spacer(minLength: 0)
                Button { model.numericFilter = nil } label: { Image(systemName: "xmark.circle.fill").frame(width: 44, height: 44) }
                    .accessibilityLabel("Clear numeric filter").accessibilityIdentifier("clearNumericFilter")
            }.font(.caption.weight(.semibold)).padding(.leading, 12)
                .background(ChemobStyle.lime.opacity(0.4), in: RoundedRectangle(cornerRadius: 12))
                .accessibilityIdentifier("activeNumericFilter")
        }
    }
    private func range(_ filter: NumericDescriptorFilter) -> String {
        func text(_ number: Double) -> String { String(number).replacingOccurrences(of: ".", with: locale.decimalSeparator ?? ".") }
        switch (filter.minimum, filter.maximum) {
        case let (.some(low), .some(high)): return low == high ? "= " + text(low) : text(low) + " – " + text(high)
        case let (.some(low), nil): return "≥ " + text(low)
        case let (nil, .some(high)): return "≤ " + text(high)
        default: return ""
        }
    }
}

private struct NumericFilterEditor: View {
    @ObservedObject var model: Workbench
    let compact: Bool
    let close: () -> Void
    @Environment(\.locale) private var locale
    @State private var descriptorID = "MolWt"
    @State private var minimum = ""
    @State private var maximum = ""
    private var draft: (value: NumericDescriptorFilter?, issue: NumericFilterIssue?) {
        do {
            let lower = try NumericDescriptorFilter.parseBound(minimum, locale: locale)
            let upper = try NumericDescriptorFilter.parseBound(maximum, locale: locale)
            return (try NumericDescriptorFilter(descriptorID: descriptorID, minimum: lower, maximum: upper), nil)
        } catch let issue as NumericFilterIssue { return (nil, issue) }
        catch { return (nil, .number) }
    }
    var body: some View {
        let parsed = draft
        VStack(alignment: .leading, spacing: 4) {
            if compact {
                propertyMenu
                bounds
                HStack { message(parsed.issue); Spacer(); actions(parsed.value) }
            } else {
                HStack(spacing: 8) {
                    propertyMenu
                    bounds.frame(width: 224)
                    actions(parsed.value)
                }
                message(parsed.issue)
            }
        }.padding(8).background(.white, in: RoundedRectangle(cornerRadius: 12))
            .onAppear {
                if let applied = model.numericFilter {
                    descriptorID = applied.descriptorID
                    let separator = locale.decimalSeparator ?? "."
                    func text(_ value: Double?) -> String {
                        value.map { String($0).replacingOccurrences(of: ".", with: separator) } ?? ""
                    }
                    minimum = text(applied.minimum); maximum = text(applied.maximum)
                }
            }
    }
    private var propertyMenu: some View {
        let catalog = Chemistry.descriptorCatalog
        let groups = catalog.reduce(into: [String]()) { if !$0.contains($1.group) { $0.append($1.group) } }
        let chosen = catalog.first { $0.id == descriptorID }
        return Menu {
            ForEach(groups, id: \.self) { group in
                Menu(LocalizedStringKey(group)) {
                    ForEach(catalog.filter { $0.group == group }) { spec in
                        Button(LocalizedStringKey(spec.title)) { descriptorID = spec.id }
                    }
                }
            }
        } label: {
            HStack(spacing: 6) {
                VStack(alignment: .leading, spacing: 2) {
                    Text(LocalizedStringKey(chosen?.title ?? "Choose a property.")).font(.subheadline.weight(.semibold)).lineLimit(1)
                    if let unit = chosen?.unit, !unit.isEmpty { Text(unit).font(.caption2).foregroundStyle(ChemobStyle.muted) }
                }
                Spacer(minLength: 0)
                Image(systemName: "chevron.down").font(.caption)
            }.frame(maxWidth: .infinity, minHeight: 44, alignment: .leading).contentShape(Rectangle())
        }.accessibilityLabel("Filter property").accessibilityIdentifier("filterDescriptor")
    }
    private var bounds: some View {
        HStack(spacing: 6) {
            TextField("Minimum", text: $minimum).accessibilityIdentifier("filterMinimum")
            Text("to").font(.caption).foregroundStyle(ChemobStyle.muted)
            TextField("Maximum", text: $maximum).accessibilityIdentifier("filterMaximum")
        }.font(.subheadline.monospacedDigit()).textFieldStyle(.roundedBorder)
            .keyboardType(.numbersAndPunctuation).textInputAutocapitalization(.never).autocorrectionDisabled()
            .frame(minHeight: 44)
    }
    private func message(_ issue: NumericFilterIssue?) -> some View {
        Text(LocalizedStringKey(issue?.rawValue ?? "Undefined values excluded"))
            .font(.caption2).foregroundStyle(issue == nil ? ChemobStyle.muted : ChemobStyle.warning)
            .lineLimit(1).frame(minHeight: 16).accessibilityIdentifier("numericFilterMessage")
    }
    private func actions(_ value: NumericDescriptorFilter?) -> some View {
        HStack(spacing: 4) {
            Button { model.numericFilter = nil; close() } label: { Image(systemName: "xmark.circle").frame(width: 44, height: 44) }
                .accessibilityLabel("Clear numeric filter").accessibilityIdentifier("clearNumericFilter")
            Button("Apply") { if let value { model.numericFilter = value; close() } }
                .font(.subheadline.bold()).padding(.horizontal, 12).frame(minHeight: 44)
                .background(ChemobStyle.lime, in: RoundedRectangle(cornerRadius: 10))
                .disabled(value == nil).accessibilityIdentifier("applyNumericFilter")
        }
    }
}

/// Fixed iPad work surface: paginated library above, selected molecule below.
struct InlineWorkbench: View {
    @ObservedObject var model: Workbench
    let add: () -> Void
    let open: () -> Void
    let export: () -> Void
    @State private var category = 0
    @State private var page = 0
    @State private var selectedID: Int?
    @State private var filtersVisible = false
    @FocusState private var searching: Bool

    private var filtered: [CompoundInput] {
        model.visible.filter { input in
            category == 0 || (category == 1 ? model.result(for: input)?.isValid == true : model.result(for: input)?.isValid == false)
        }
    }
    var body: some View {
        GeometryReader { geometry in
            let values = filtered
            let topHeight = geometry.size.height * 0.5
            let capacity = max(1, Int((topHeight - 232 - (filtersVisible ? 96 : model.numericFilter == nil ? 0 : 44)) / 62))
            let pages = max(1, (values.count + capacity - 1) / capacity)
            let current = min(page, pages - 1)
            let start = current * capacity
            let shown = Array(values.dropFirst(start).prefix(capacity))
            let active = values.first { $0.id == selectedID } ?? shown.first
            VStack(spacing: 12) {
                VStack(spacing: 8) {
                    controls
                    if filtersVisible { NumericFilterEditor(model: model, compact: false) { filtersVisible = false } }
                    else if model.numericFilter != nil { NumericFilterSummary(model: model) }
                    Picker("Category", selection: $category) {
                        Text("All").tag(0)
                        Text("Valid").tag(1)
                        Text("Needs review").tag(2)
                    }.pickerStyle(.segmented).frame(minHeight: 44).accessibilityIdentifier("compoundCategory")
                    HStack {
                        Text("COMPOUND").frame(maxWidth: .infinity, alignment: .leading)
                        Text("MW").frame(width: 72, alignment: .trailing)
                        Text("cLogP").frame(width: 76, alignment: .trailing)
                        Text("TPSA / Å²").frame(width: 90, alignment: .trailing)
                    }.font(.caption2.bold()).foregroundStyle(ChemobStyle.muted).padding(.horizontal, 16)
                    VStack(spacing: 0) {
                        ForEach(shown) { input in
                            Button { selectedID = input.id; searching = false } label: {
                                InlineCompoundRow(input: input, result: model.result(for: input), selected: active?.id == input.id)
                            }.buttonStyle(CompoundRowStyle()).accessibilityIdentifier("compound-\(input.id)")
                                .accessibilityAddTraits(active?.id == input.id ? .isSelected : [])
                        }
                        if shown.isEmpty {
                            Text("No matching compounds").foregroundStyle(ChemobStyle.muted).frame(maxWidth: .infinity, maxHeight: .infinity)
                        }
                    }.background(.white, in: RoundedRectangle(cornerRadius: 12))
                    Spacer(minLength: 0)
                    HStack {
                        Group {
                            if values.count == 1 { Text("1 match") }
                            else { Text("\(values.count) matches") }
                        }.font(.caption).foregroundStyle(ChemobStyle.muted)
                        Spacer()
                        Button { page = current - 1; selectedID = nil; searching = false } label: {
                            Image(systemName: "chevron.left").frame(width: 44, height: 44)
                        }.disabled(current == 0).accessibilityLabel("Previous page")
                        Text("\(current + 1) / \(pages)").font(.caption.monospacedDigit())
                        Button { page = current + 1; selectedID = nil; searching = false } label: {
                            Image(systemName: "chevron.right").frame(width: 44, height: 44)
                        }.disabled(current + 1 >= pages).accessibilityLabel("Next page")
                    }
                }.frame(height: max(0, topHeight - 12))
                if let active {
                    InlineCompoundDetail(input: active, result: model.result(for: active), model: model)
                        .frame(maxWidth: .infinity, maxHeight: .infinity)
                } else {
                    VStack(spacing: 12) {
                        Image(systemName: "magnifyingglass").font(.title)
                        Text("No matching compounds").font(.headline)
                    }.foregroundStyle(ChemobStyle.muted).frame(maxWidth: .infinity, maxHeight: .infinity).background(.white, in: RoundedRectangle(cornerRadius: 20))
                }
            }
        }.onChange(of: model.query) { page = 0; selectedID = nil }
            .onChange(of: category) { page = 0; selectedID = nil }
            .onChange(of: model.numericFilter) { page = 0; selectedID = nil; searching = false }
    }
    private var controls: some View {
        VStack(spacing: 8) {
            HStack(spacing: 12) {
                Text("Chemob").font(.system(.title2, design: .rounded).bold())
                Text(model.inputs.count, format: .number).font(.headline.monospacedDigit()).foregroundStyle(ChemobStyle.muted)
                Spacer()
                Button(action: add) { Label("Add SMILES", systemImage: "plus").padding(.horizontal, 14).frame(minHeight: 44) }
                    .buttonStyle(WorkbenchButtonStyle(primary: true)).accessibilityIdentifier("pasteSMILES")
                Button(action: open) { Label("Import", systemImage: "folder").padding(.horizontal, 12).frame(minHeight: 44) }
                    .buttonStyle(WorkbenchButtonStyle(primary: false)).accessibilityIdentifier("importFile")
                if model.busy { Button("Cancel") { model.cancel() }.frame(minHeight: 44) }
                else {
                    CSVExportControl(model: model, action: export, iconOnly: true)
                }
            }.font(.subheadline.bold())
            HStack(spacing: 8) {
                Image(systemName: "magnifyingglass").accessibilityHidden(true)
                TextField("Search name or SMILES", text: $model.query).focused($searching)
                    .textInputAutocapitalization(.never).autocorrectionDisabled().submitLabel(.search).onSubmit { searching = false }
                    .accessibilityIdentifier("compoundSearch")
                if !model.query.isEmpty {
                    Button { model.query = "" } label: { Image(systemName: "xmark.circle.fill").frame(width: 44, height: 44) }.accessibilityLabel("Clear search")
                }
                NumericFilterToggle(active: model.numericFilter != nil, expanded: $filtersVisible)
            }.padding(.horizontal, 14).frame(minHeight: 44).background(.white, in: RoundedRectangle(cornerRadius: 12))
        }
    }
}

private struct InlineCompoundRow: View {
    let input: CompoundInput
    let result: AnalysisRow?
    let selected: Bool
    @Environment(\.locale) private var locale
    var body: some View {
        HStack(spacing: 16) {
            VStack(alignment: .leading, spacing: 3) {
                Text(input.name).font(.subheadline.weight(.semibold)).lineLimit(1)
                Text(result?.formula ?? input.smiles).font(.caption.monospaced()).foregroundStyle(ChemobStyle.muted).lineLimit(1)
            }.frame(maxWidth: .infinity, alignment: .leading)
            number(result?.mw, width: 72)
            number(result?.logp, width: 76)
            number(result?.tpsa, width: 90)
        }.padding(.horizontal, 16).frame(height: 62)
            .background(selected ? ChemobStyle.lime.opacity(0.3) : Color.clear)
            .overlay(alignment: .leading) { if selected { ChemobStyle.ink.frame(width: 3).padding(.vertical, 10) } }
            .contentShape(Rectangle())
    }
    private func number(_ value: Double?, width: CGFloat) -> some View {
        Text(DescriptorFormatting.text(value, digits: 2, locale: locale)).font(.subheadline.monospacedDigit()).frame(width: width, alignment: .trailing)
    }
}

private struct InlineCompoundDetail: View {
    let input: CompoundInput
    let result: AnalysisRow?
    @ObservedObject var model: Workbench
    @State private var graph: ChemicalDepiction?
    @State private var graphSource = ""
    @State private var loadingGraph = false
    var body: some View {
        GeometryReader { geometry in
            VStack(alignment: .leading, spacing: 8) {
                HStack(alignment: .top) {
                    VStack(alignment: .leading, spacing: 4) {
                        Text(input.name).font(.system(.title2, design: .rounded).bold()).lineLimit(1).accessibilityIdentifier("selectedCompoundName")
                        Text(result?.formula ?? "SMILES").font(.subheadline.monospaced()).foregroundStyle(ChemobStyle.muted)
                    }
                    Spacer()
                    Label(LocalizedStringKey(result == nil ? "Not calculated" : (result?.isValid == true ? "Ready" : "Needs review")), systemImage: result == nil ? "hourglass" : (result?.isValid == true ? "checkmark.circle" : "exclamationmark.triangle"))
                        .font(.caption.bold()).padding(10).background(result?.isValid == false ? ChemobStyle.warning.opacity(0.1) : ChemobStyle.lime, in: Capsule())
                        .accessibilityIdentifier("selectedCompoundStatus")
                }
                if let result, result.isValid {
                    HStack(alignment: .top, spacing: 20) {
                        VStack(spacing: 8) {
                            if let graph, graphSource == input.smiles { ChemicalStructurePreview(drawing: graph) }
                            else if loadingGraph { ProgressView().frame(maxWidth: .infinity, maxHeight: .infinity) }
                            else {
                                Label("Preview unavailable", systemImage: "square.dashed")
                                    .font(.caption).foregroundStyle(ChemobStyle.muted).multilineTextAlignment(.center)
                                    .frame(maxWidth: .infinity, maxHeight: .infinity)
                            }
                            Text(LocalizedStringKey(graph?.stereochemistryOmitted == true ? "2D structure · some stereochemistry not shown" : graph == nil ? "Preview unavailable" : "2D structure")).font(.caption2).foregroundStyle(ChemobStyle.muted).multilineTextAlignment(.center)
                        }.frame(width: geometry.size.width * 0.27)
                        descriptorGrid(result, columns: geometry.size.width >= 1000 ? 4 : 3)
                    }.frame(maxHeight: .infinity)
                } else {
                    Label(LocalizedStringKey(result == nil ? "This molecule has not been calculated yet." : "Invalid or unsupported structure. The original SMILES is preserved."), systemImage: result == nil ? "hourglass" : "exclamationmark.triangle")
                        .foregroundStyle(result == nil ? ChemobStyle.muted : ChemobStyle.warning).frame(maxWidth: .infinity, maxHeight: .infinity, alignment: .leading)
                }
                HStack(alignment: .top, spacing: 12) {
                    Text("SMILES").font(.caption.bold()).foregroundStyle(ChemobStyle.muted)
                    Text(input.smiles).font(.caption.monospaced()).lineLimit(2).textSelection(.enabled).frame(maxWidth: .infinity, alignment: .leading)
                }
                Divider()
                HStack {
                    Label("On device", systemImage: "checkmark.shield")
                    Spacer()
                    if model.busy { Text("Processed \(model.results.count) of \(model.inputs.count)") }
                    else if model.exporting { Text("Exporting CSV") }
                    else if model.reusedCount > 0 { Text("Reused \(model.reusedCount) results") }
                    else { Text("Calculated locally") }
                    if let elapsed = model.elapsed { Text(elapsed, format: .number.precision(.fractionLength(2))) + Text(" s") }
                }.font(.caption).foregroundStyle(ChemobStyle.muted).accessibilityIdentifier("calculationStatus")
            }.padding(16).frame(maxWidth: .infinity, maxHeight: .infinity, alignment: .topLeading)
                .background(.white, in: RoundedRectangle(cornerRadius: 20))
                .accessibilityElement(children: .contain).accessibilityIdentifier("inlineCompoundDetail")
        }.task(id: input.smiles + (result?.status ?? "pending")) {
            graph = nil
            loadingGraph = result?.isValid == true
            guard result?.isValid == true else { return }
            let value = await ChemistryStore.shared.drawing(input.smiles)
            if !Task.isCancelled { graphSource = input.smiles; graph = value; loadingGraph = false }
        }
    }
    private func descriptorGrid(_ row: AnalysisRow, columns: Int) -> some View {
        DescriptorBrowser(row: row, columns: columns).frame(maxWidth: .infinity, alignment: .topLeading)
    }
}

private struct DescriptorBrowser: View {
    let row: AnalysisRow
    let columns: Int
    @State private var group = "Essentials"
    @State private var page = 0
    @State private var informationPresented = false
    @Environment(\.locale) private var locale
    private var legacyValues: [Double?] {
        [row.mw, row.exact_mass, row.logp, row.tpsa, row.hbd.map(Double.init), row.hba.map(Double.init), row.heavy_atoms.map(Double.init), row.bonds.map(Double.init), row.formal_charge.map(Double.init), row.fragments.map(Double.init), row.rotatable_bonds.map(Double.init), row.fraction_csp3]
    }
    var body: some View {
        let catalog = Chemistry.descriptorCatalog
        let groups = catalog.reduce(into: [String]()) { if !$0.contains($1.group) { $0.append($1.group) } }
        let indices = catalog.indices.filter { catalog[$0].group == group }
        let capacity = columns >= 4 ? 8 : 12
        let pages = max(1, (indices.count + capacity - 1) / capacity)
        let current = min(page, pages - 1)
        let shown = Array(indices.dropFirst(current * capacity).prefix(capacity))
        let values = row.descriptor_values ?? legacyValues
        VStack(alignment: .leading, spacing: 8) {
            HStack(spacing: 4) {
                Picker("Properties", selection: $group) {
                    ForEach(groups, id: \.self) { Text(LocalizedStringKey($0)).tag($0) }
                }.pickerStyle(.menu).font(.caption.bold()).frame(minWidth: 150, alignment: .leading).accessibilityIdentifier("descriptorGroup")
                Spacer(minLength: 0)
                if shown.contains(where: { ["Ipc", "AvgIpc"].contains(catalog[$0].id) }) {
                    Button { informationPresented = true } label: {
                        Image(systemName: "info.circle").frame(width: 44, height: 44)
                    }.accessibilityLabel("About graph information")
                        .accessibilityIdentifier("graphInformationInfo")
                        .popover(isPresented: $informationPresented) {
                            Text("Stable information compatibility")
                                .font(.callout).padding(20).frame(idealWidth: 320)
                                .fixedSize(horizontal: false, vertical: true)
                                .accessibilityIdentifier("graphInformationExplanation")
                                .presentationCompactAdaptation(.popover)
                        }
                }
                if pages > 1 {
                    Button { page = current - 1 } label: { Image(systemName: "chevron.left").frame(width: 44, height: 44) }
                        .disabled(current == 0).accessibilityLabel("Previous properties")
                    Text("\(current + 1) / \(pages)").font(.caption.monospacedDigit())
                    Button { page = current + 1 } label: { Image(systemName: "chevron.right").frame(width: 44, height: 44) }
                        .disabled(current + 1 >= pages).accessibilityLabel("Next properties")
                }
            }.frame(minHeight: 44)
            LazyVGrid(columns: Array(repeating: GridItem(.flexible(), alignment: .leading), count: columns), alignment: .leading, spacing: 8) {
                ForEach(shown, id: \.self) { index in
                    let spec = catalog[index]
                    let value = values.indices.contains(index) ? values[index] : nil
                    VStack(alignment: .leading, spacing: 4) {
                        Text(LocalizedStringKey(spec.title)).font(.caption2).foregroundStyle(ChemobStyle.muted).lineLimit(2)
                        HStack(alignment: .firstTextBaseline, spacing: 4) {
                            Text(DescriptorFormatting.text(value, digits: spec.digits, locale: locale))
                                .font(.system(.title3, design: .rounded).weight(.semibold)).monospacedDigit().lineLimit(1).minimumScaleFactor(0.7)
                                .accessibilityIdentifier("descriptor-" + spec.id)
                            if !spec.unit.isEmpty { Text(spec.unit).font(.caption2).foregroundStyle(ChemobStyle.muted).lineLimit(1) }
                        }
                    }.frame(maxWidth: .infinity, minHeight: 52, alignment: .topLeading)
                }
            }
        }.onChange(of: group) { page = 0 }
    }
}
