import { useState } from 'react';
import { motion } from 'framer-motion';

const ContactSection = () => {
  const [formData, setFormData] = useState({
    name: '',
    email: '',
    institution: '',
    message: ''
  });
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [isSubmitted, setIsSubmitted] = useState(false);

  const handleSubmit = async (e) => {
    e.preventDefault();
    setIsSubmitting(true);
    
    // Simulate form submission
    await new Promise(resolve => setTimeout(resolve, 1500));
    
    setIsSubmitting(false);
    setIsSubmitted(true);
  };

  const handleChange = (e) => {
    setFormData(prev => ({
      ...prev,
      [e.target.name]: e.target.value
    }));
  };

  return (
    <section id="contact" className="relative w-full min-h-screen bg-[#0A0A0A] py-20 px-4 sm:px-8 lg:px-16">
      {/* Magenta Bloom Effect */}
      <div 
        className="absolute right-[-10vw] bottom-[-10vh] w-[50vw] h-[50vh] pointer-events-none"
        style={{
          background: 'radial-gradient(ellipse at center, rgba(255, 45, 141, 0.1) 0%, transparent 60%)',
          filter: 'blur(80px)'
        }}
      />

      <div className="relative z-10 max-w-7xl mx-auto">
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-12 lg:gap-20">
          {/* Left Column - Copy */}
          <motion.div
            initial={{ opacity: 0, y: 40 }}
            whileInView={{ opacity: 1, y: 0 }}
            viewport={{ once: true, margin: "-100px" }}
            transition={{ duration: 0.6 }}
          >
            <h2 className="font-display text-[clamp(32px,4vw,52px)] font-bold text-[#F2F4F8] leading-tight">
              Ready to analyze
              <span className="block text-transparent bg-clip-text bg-gradient-to-r from-[#00F0FF] to-[#FF2D8D]">
                your network?
              </span>
            </h2>
            
            <p className="mt-6 text-[#A6AEB8] text-lg leading-relaxed max-w-md">
              Send us your dataset details. We'll reply with a pipeline plan, 
              estimated runtime, and analysis recommendations.
            </p>

            {/* Features */}
            <div className="mt-10 space-y-4">
              {[
                'Free initial consultation',
                'Custom pipeline configuration',
                'Detailed analysis report',
                'Post-analysis support'
              ].map((feature, i) => (
                <motion.div
                  key={feature}
                  className="flex items-center gap-3"
                  initial={{ opacity: 0, x: -20 }}
                  whileInView={{ opacity: 1, x: 0 }}
                  viewport={{ once: true }}
                  transition={{ duration: 0.4, delay: 0.3 + i * 0.1 }}
                >
                  <div className="w-5 h-5 rounded-full bg-[#00F0FF]/20 flex items-center justify-center">
                    <svg className="w-3 h-3 text-[#00F0FF]" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={3} d="M5 13l4 4L19 7" />
                    </svg>
                  </div>
                  <span className="text-[#F2F4F8]">{feature}</span>
                </motion.div>
              ))}
            </div>

            {/* Contact Info */}
            <div className="mt-12 pt-8 border-t border-white/10">
              <div className="flex items-center gap-4">
                <div className="w-12 h-12 rounded-xl bg-[#00F0FF]/10 flex items-center justify-center">
                  <svg className="w-6 h-6 text-[#00F0FF]" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} 
                          d="M3 8l7.89 5.26a2 2 0 002.22 0L21 8M5 19h14a2 2 0 002-2V7a2 2 0 00-2-2H5a2 2 0 00-2 2v10a2 2 0 002 2z" />
                  </svg>
                </div>
                <div>
                  <div className="text-sm text-[#A6AEB8]">Email us at</div>
                  <div className="font-mono-data text-[#F2F4F8]">contact@dice-duan.org</div>
                </div>
              </div>
            </div>
          </motion.div>

          {/* Right Column - Form */}
          <motion.div
            initial={{ opacity: 0, x: 40 }}
            whileInView={{ opacity: 1, x: 0 }}
            viewport={{ once: true, margin: "-100px" }}
            transition={{ duration: 0.6, delay: 0.1 }}
          >
            <div className="glass-card p-8">
              {isSubmitted ? (
                <motion.div
                  className="text-center py-12"
                  initial={{ opacity: 0, scale: 0.95 }}
                  animate={{ opacity: 1, scale: 1 }}
                  transition={{ duration: 0.4 }}
                >
                  <div className="w-20 h-20 mx-auto mb-6 rounded-full bg-[#00F0FF]/10 flex items-center justify-center">
                    <svg className="w-10 h-10 text-[#00F0FF]" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
                    </svg>
                  </div>
                  <h3 className="font-display text-2xl font-bold text-[#F2F4F8] mb-2">
                    Request received!
                  </h3>
                  <p className="text-[#A6AEB8]">
                    We'll get back to you within 24 hours.
                  </p>
                  <button
                    onClick={() => {
                      setIsSubmitted(false);
                      setFormData({ name: '', email: '', institution: '', message: '' });
                    }}
                    className="mt-6 text-[#00F0FF] hover:underline"
                  >
                    Send another request
                  </button>
                </motion.div>
              ) : (
                <form onSubmit={handleSubmit} className="space-y-5">
                  <div>
                    <label className="block text-sm text-[#A6AEB8] mb-2">Name</label>
                    <input
                      type="text"
                      name="name"
                      value={formData.name}
                      onChange={handleChange}
                      required
                      className="w-full px-4 py-3 rounded-xl bg-white/[0.03] border border-white/[0.08]
                                 text-[#F2F4F8] placeholder-[#A6AEB8]/50
                                 focus:outline-none focus:border-[#00F0FF]/50 focus:ring-1 focus:ring-[#00F0FF]/30
                                 transition-all duration-200"
                      placeholder="Your name"
                    />
                  </div>

                  <div>
                    <label className="block text-sm text-[#A6AEB8] mb-2">Email</label>
                    <input
                      type="email"
                      name="email"
                      value={formData.email}
                      onChange={handleChange}
                      required
                      className="w-full px-4 py-3 rounded-xl bg-white/[0.03] border border-white/[0.08]
                                 text-[#F2F4F8] placeholder-[#A6AEB8]/50
                                 focus:outline-none focus:border-[#00F0FF]/50 focus:ring-1 focus:ring-[#00F0FF]/30
                                 transition-all duration-200"
                      placeholder="you@institution.edu"
                    />
                  </div>

                  <div>
                    <label className="block text-sm text-[#A6AEB8] mb-2">Institution</label>
                    <input
                      type="text"
                      name="institution"
                      value={formData.institution}
                      onChange={handleChange}
                      className="w-full px-4 py-3 rounded-xl bg-white/[0.03] border border-white/[0.08]
                                 text-[#F2F4F8] placeholder-[#A6AEB8]/50
                                 focus:outline-none focus:border-[#00F0FF]/50 focus:ring-1 focus:ring-[#00F0FF]/30
                                 transition-all duration-200"
                      placeholder="University or company"
                    />
                  </div>

                  <div>
                    <label className="block text-sm text-[#A6AEB8] mb-2">Message</label>
                    <textarea
                      name="message"
                      value={formData.message}
                      onChange={handleChange}
                      rows={4}
                      className="w-full px-4 py-3 rounded-xl bg-white/[0.03] border border-white/[0.08]
                                 text-[#F2F4F8] placeholder-[#A6AEB8]/50 resize-none
                                 focus:outline-none focus:border-[#00F0FF]/50 focus:ring-1 focus:ring-[#00F0FF]/30
                                 transition-all duration-200"
                      placeholder="Tell us about your dataset and research goals..."
                    />
                  </div>

                  <button
                    type="submit"
                    disabled={isSubmitting}
                    className="w-full py-4 px-6 rounded-xl bg-[#00F0FF]/10 border border-[#00F0FF]/40
                               text-[#00F0FF] font-medium
                               hover:bg-[#00F0FF]/20 hover:border-[#00F0FF]/60
                               disabled:opacity-50 disabled:cursor-not-allowed
                               transition-all duration-200 flex items-center justify-center gap-2"
                  >
                    {isSubmitting ? (
                      <>
                        <svg className="w-5 h-5 animate-spin" fill="none" viewBox="0 0 24 24">
                          <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                          <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z" />
                        </svg>
                        Sending...
                      </>
                    ) : (
                      <>
                        Request analysis
                        <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M17 8l4 4m0 0l-4 4m4-4H3" />
                        </svg>
                      </>
                    )}
                  </button>
                </form>
              )}
            </div>
          </motion.div>
        </div>

        {/* Footer */}
        <motion.footer
          className="mt-20 pt-8 border-t border-white/10"
          initial={{ opacity: 0 }}
          whileInView={{ opacity: 1 }}
          viewport={{ once: true }}
          transition={{ duration: 0.5, delay: 0.3 }}
        >
          <div className="flex flex-col md:flex-row items-center justify-between gap-4">
            <div className="flex items-center gap-2">
              <span className="font-display font-bold text-[#F2F4F8]">DiCE-Duan</span>
              <span className="text-[#A6AEB8]">—</span>
              <span className="text-sm text-[#A6AEB8]">Built for high-resolution network biology</span>
            </div>
            
            <div className="flex items-center gap-6">
              <a href="#" className="text-sm text-[#A6AEB8] hover:text-[#00F0FF] transition-colors">
                Documentation
              </a>
              <a href="#" className="text-sm text-[#A6AEB8] hover:text-[#00F0FF] transition-colors">
                GitHub
              </a>
              <a href="#" className="text-sm text-[#A6AEB8] hover:text-[#00F0FF] transition-colors">
                Citation
              </a>
            </div>
          </div>
          
          <div className="mt-4 text-center text-xs text-[#A6AEB8]/50">
            © {new Date().getFullYear()} DiCE-Duan Project. All rights reserved.
          </div>
        </motion.footer>
      </div>
    </section>
  );
};

export default ContactSection;
