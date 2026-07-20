import sys
from pathlib import Path
import step5_db_search as db

DB_PATH = "viral9_refs.jsonl"


def get_payload(path_or_payload):
    p = Path(path_or_payload)

    # If user gives payload text file
    if p.exists() and p.suffix.lower() == ".txt":
        return p.read_text(encoding="utf-8").strip()

    # If user gives QR image
    if p.exists() and p.suffix.lower() in [".png", ".jpg", ".jpeg"]:
        # Try pyzbar first
        try:
            from PIL import Image
            from pyzbar.pyzbar import decode

            img = Image.open(str(p))
            decoded = decode(img)

            if decoded:
                return decoded[0].data.decode("utf-8").strip()
        except Exception as e:
            print(f"pyzbar failed: {e}")

        # Fallback to OpenCV
        try:
            import cv2

            img = cv2.imread(str(p))
            detector = cv2.QRCodeDetector()
            payload, points, _ = detector.detectAndDecode(img)

            if payload:
                return payload.strip()
        except Exception as e:
            print(f"OpenCV failed: {e}")

        raise ValueError("Could not decode QR image using pyzbar or OpenCV.")

    # If user directly pastes VGQR1 payload
    return path_or_payload.strip()


def main():
    if len(sys.argv) < 2:
        print("Usage:")
        print("  py qr_decode_search.py qr_image.png")
        print("  py qr_decode_search.py payload.txt")
        print("  py qr_decode_search.py VGQR1....")
        return

    query_input = sys.argv[1]
    payload = get_payload(query_input)

    query = db.decode_payload(payload)
    refs = db.load_db_jsonl(DB_PATH)

    results = []
    for ref in refs:
        mismatch = db.check_compat(query, ref)
        if mismatch:
            raise ValueError(f"Query/reference parameter mismatch: {mismatch}")

        sim = db.minhash_jaccard_est(query["sketch"], ref["sketch"])
        results.append((ref["name"], sim))

    results.sort(key=lambda x: x[1], reverse=True)

    top_name, top_score = results[0]
    second_name, second_score = results[1]
    gap = top_score - second_score
    confidence = db.confidence_label(top_score, gap)

    print("\n=== QR Decode + Viral DB Search Result ===")
    print(f"Top match: {top_name}")
    print(f"Similarity: {top_score:.4f}")
    print(f"Second match: {second_name}")
    print(f"Second similarity: {second_score:.4f}")
    print(f"Confidence gap: {gap:.4f}")
    print(f"Confidence: {confidence}")

    print("\nTop 5 matches:")
    for i, (name, score) in enumerate(results[:5], start=1):
        print(f"{i}. {name}: {score:.4f}")


if __name__ == "__main__":
    main()
