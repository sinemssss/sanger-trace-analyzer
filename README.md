# Sanger Trace Analyzer 🧬

Genetik tanı laboratuvarlarında Sanger dizileme (`.ab1`) sonuçlarının
**manuel** pik-pik incelenmesini hızlandırmak için geliştirilmiş bir Python
aracı. Her pozisyondaki 4 kanalı (A/C/G/T) karşılaştırıp, tek pik (normal)
ile çift pik (heterozigot / mixed / mutasyon şüphesi) olan pozisyonları
otomatik ayırt eder — biyolog yalnızca işaretlenen pozisyonlara bakar.

![örnek çıktı](examples/kromatogram_zoom.png)

## Neden

Bir genetik tanı merkezinde staj sırasında, YMİK ve benzeri analizlerde her
pikin manuel kontrol edildiğini gözlemledim. Bu araç, günün sonunda çıkan çok
sayıda `.ab1` dosyasını otomatik tarayarak sadece "dikkat gerektiren"
pozisyonları öne çıkarır. STR, kimerizm ve kanser mutasyon taramaları gibi
heterozigot/mixed pik tespitinin önemli olduğu senaryolara uyarlanabilir.

**Not:** Bu bir portföy/prototip projesidir; klinik kullanım için gerçek
laboratuvar verisiyle kapsamlı doğrulama ve kalibrasyon gerekir (bkz.
[Sınırlamalar](#önemli-notlar--sınırlamalar)).

## Öne çıkan özellik: bağımsız `.ab1` okuyucu

`.ab1` (ABIF) dosyaları genellikle Biopython ile okunur. Bu projede,
formatı **sıfırdan**, kamuya açık ABIF spesifikasyonuna göre okuyan bağımsız
bir parser yazdım (`abif_reader.py`) — harici bioinformatik kütüphanesine
bağımlılık yok, sadece `pandas` ve `matplotlib`.

## Kurulum

```bash
git clone https://github.com/<kullanici-adin>/sanger-trace-analyzer.git
cd sanger-trace-analyzer
pip install -r requirements.txt
# veya CLI komutu olarak kurmak için:
pip install -e .
```

## Kullanım

```bash
# Doğrudan modül olarak
python3 -m sanger_trace_analyzer.main examples/test_sample.ab1 --plot kromatogram.png

# pip install -e . sonrası CLI komutu olarak
sanger-analyze examples/test_sample.ab1 --plot kromatogram.png

# Referans (yabani tip) dizisiyle karşılaştırma
sanger-analyze ornek.ab1 --ref ACGTACGT... --threshold 0.20
```

Çıktı:
- Terminalde özet (kaç pozisyon normal / şüpheli / düşük sinyal)
- `rapor.csv`: her pozisyon için pik yükseklikleri ve oranları
- (opsiyonel) kromatogram görseli + şüpheli pozisyonların yakınlaştırılmış görünümü

## Nasıl çalışır

1. **`abif_reader.py`** — `.ab1` (ABIF) ikili dosyasını sıfırdan okur:
   ham trace kanalları (DATA9-12), baz çağrıları (PBAS) ve pik konumları
   (PLOC) çıkarılır.
2. **`sanger_analyzer.py`** — her baz pozisyonunda 4 kanalın o noktadaki
   pik yüksekliğine bakar; en yüksek iki kanal "primary" / "secondary"
   olarak alınır. `secondary/primary` oranı bir eşiği (varsayılan **0.25**)
   aşarsa ve gürültü tabanının belirgin üzerindeyse pozisyon
   `mixed_candidate` olarak işaretlenir. Gürültü tabanı, ortalama/std yerine
   **medyan + MAD** (median absolute deviation) ile hesaplanır — bu, pik
   yüksekliklerinden etkilenmeyen daha sağlam bir yöntemdir.
3. **`plot_trace.py`** — tam kromatogramı ve şüpheli pozisyonların
   yakınlaştırılmış görünümünü çizer.
4. **`main.py`** — komut satırı arayüzü.

## Test

Gerçek laboratuvar verisi olmadan da doğrulanabilmesi için,
`examples/make_fake_ab1.py` bilinçli olarak yerleştirilmiş çift-pik
pozisyonları içeren **sentetik** bir `.ab1` dosyası üretir (gerçek biyolojik
veri değildir). `tests/test_analyzer.py` bu sentetik veriyle algoritmanın
doğruluğunu kontrol eder:

```bash
python3 tests/test_analyzer.py -v
# ya da pytest kuruluysa:
pytest tests/ -v
```

## Proje yapısı

```
sanger-trace-analyzer/
├── sanger_trace_analyzer/   # ana paket
│   ├── abif_reader.py       # .ab1 dosya okuyucu (Biopython gerektirmez)
│   ├── sanger_analyzer.py   # pik analiz algoritması
│   ├── plot_trace.py        # görselleştirme
│   └── main.py               # CLI
├── examples/
│   ├── make_fake_ab1.py     # test için sentetik .ab1 üretici
│   ├── test_sample.ab1      # sentetik örnek dosya
│   └── kromatogram*.png     # örnek çıktılar
├── tests/
│   └── test_analyzer.py
├── requirements.txt
├── setup.py
└── README.md
```

## Önemli notlar / sınırlamalar

- **Eşik değerleri (`--threshold`, `noise_floor_k`) kalibrasyon gerektirir.**
  Bu sürümdeki varsayılanlar (0.25 oran eşiği) literatürde sık kullanılan
  başlangıç değerleridir, ancak cihaz/kit/panel'e göre gerçek `.ab1`
  dosyalarıyla ayarlanmalıdır (yanlış pozitif/negatif oranını görmek için
  bilinen sonuçlu birkaç örnekle test edin).
- Şu an tek nokta mutasyonlarına / heterozigot pozisyonlara odaklanıyor;
  insersiyon/delesyon (indel) sonrası kayan çift pik desenleri için ayrı
  bir mantık eklenmesi gerekir — doğal bir sonraki adım.
- Referans karşılaştırması şu anda basit indeks-bazlı; indel içeren
  karşılaştırmalar için hizalama (alignment) eklenmeli.
- **Gerçek hasta/örnek verisi bu repoya asla commit edilmemelidir**
  (`.gitignore` `*.ab1` dosyalarını hariç tutar, sadece sentetik
  `examples/test_sample.ab1` istisnadır).

## Sonraki adımlar (öneriler)

- Gerçek `.ab1` dosyalarıyla eşik kalibrasyonu (ROC eğrisi ile optimum eşik)
- İndel-farkında hizalama (örn. basit Needleman-Wunsch) ile referans karşılaştırması
- Toplu (batch) mod: bir klasördeki tüm `.ab1` dosyalarını tarayıp özet rapor üretme
- STR analizi için: tekrar sayısı tahmini modülü
- Kimerizm için: donör/alıcı oranı hesaplama modülü (bilinen polimorfik pozisyonlarda)

## Lisans

MIT — bkz. [LICENSE](LICENSE)
