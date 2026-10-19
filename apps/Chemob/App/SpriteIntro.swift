import SwiftUI
import ImageIO

struct SpriteIntro: View {
    @Environment(\.accessibilityReduceMotion) private var reduceMotion
    @Environment(\.scenePhase) private var scenePhase
    @State private var started = Date()
    @State private var finished = false
    private static let asset: (frames: [CGImage], duration: Double) = {
        guard let manifestURL = Bundle.main.url(forResource: "chemob-intro", withExtension: "json"),
              let data = try? Data(contentsOf: manifestURL),
              let object = try? JSONSerialization.jsonObject(with: data) as? [String: Any],
              object["schema"] as? String == "backbeat-sprite-atlas-v1",
              let width = object["frameWidth"] as? Int, let height = object["frameHeight"] as? Int,
              let origins = object["frames"] as? [[String: Any]],
              let imageURL = Bundle.main.url(forResource: "chemob-intro", withExtension: "png"),
              let source = CGImageSourceCreateWithURL(imageURL as CFURL, nil),
              let image = CGImageSourceCreateImageAtIndex(source, 0, nil),
              let duration = object["duration"] as? Double, duration.isFinite, duration > 0 else { return ([], 2.4) }
        let frames = origins.compactMap { origin -> CGImage? in
            guard let x = origin["x"] as? Int, let y = origin["y"] as? Int else { return nil }
            return image.cropping(to: CGRect(x: x, y: y, width: width, height: height))
        }
        return (frames, duration)
    }()
    var body: some View {
        TimelineView(.animation(minimumInterval: 0.1, paused: finished || reduceMotion || scenePhase != .active)) { timeline in
            let frames = Self.asset.frames
            if !frames.isEmpty {
                let progress = reduceMotion ? 1 : max(0, min(1, timeline.date.timeIntervalSince(started) / Self.asset.duration))
                let index = min(frames.count - 1, Int(progress * Double(frames.count)))
                Image(decorative: frames[index], scale: 1).resizable().interpolation(.high).scaledToFit()
            } else {
                Image(systemName: "hexagon.fill").font(.system(size: 130)).foregroundStyle(ChemobStyle.lime)
            }
        }.accessibilityHidden(true).task {
            started = Date()
            do { try await Task.sleep(for: .seconds(Self.asset.duration)); finished = true } catch { }
        }
    }
}
