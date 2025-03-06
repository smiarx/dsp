#include "../../AllPass.h"
#include "../../Buffer.h"
#include "../../Delay.h"
#include "../../Filter.h"

class Reverb
{
  public:
    static constexpr auto MaxBlockSize = 128;
    static constexpr auto a            = 0.75;
    static constexpr auto fb           = 0.57;
    static constexpr auto lpfreq = 0.2f;

  private:
    static constexpr dsp::AllPass<1> allpass1{{a}};
    static constexpr dsp::DelayLine<150> apd1{};
    static constexpr dsp::DelayLine<213, nextTo(apd1)> apd2{};
    static constexpr dsp::DelayLine<319, nextTo(apd2)> apd3{};
    static constexpr dsp::DelayLine<526, nextTo(apd3)> apd4{};

    dsp::IIRFilter<2> lpfilter {dsp::IIRFilter<2>::newButterworthLP(dsp::Signal<1>{lpfreq})};
    decltype(lpfilter)::DL lpfilterd{};

    static constexpr dsp::AllPass<2, dsp::TapFix<2182, 2525>> ap2n1{{-a, -a}};
    static constexpr dsp::DelayLine<2525> dapd1{};

    static constexpr dsp::AllPass<2, dsp::TapFix<2690, 2197>> ap2n2{{-a, -a}};
    static constexpr dsp::DelayLine<2690, nextTo(dapd1)> dapd2{};

    static constexpr dsp::AllPass<2, dsp::TapFix<2673, 2847>> ap2n3{{-a, -a}};
    static constexpr dsp::DelayLine<2847, nextTo(dapd2)> dapd3{};

    static constexpr dsp::AllPass<2, dsp::TapFix<3430, 3210>> ap2n4{{-a, -a}};
    static constexpr dsp::DelayLine<3430, nextTo(dapd3)> dapd4{};

    static constexpr dsp::DelayLine<6261, nextTo(dapd4)> feedbackd{};

  public:
    dsp::Buffer<dsp::Signal<1>, nextTo(apd4) + MaxBlockSize> buffer1;
    dsp::Buffer<dsp::Signal<2>, nextTo(feedbackd) + MaxBlockSize> buffer2;
    static constexpr auto Buffer1Size = decltype(buffer1)::Size;
    static constexpr auto Buffer2Size = decltype(buffer2)::Size;

    dsp::Signal<2> *xblock;

    void process(float **in, float **out, int count);
};
