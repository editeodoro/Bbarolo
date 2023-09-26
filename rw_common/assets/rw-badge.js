const RapidWeaverBadgeTemplate = document.createElement('template');
RapidWeaverBadgeTemplate.innerHTML = `
<style>
  [data-theme="dark"] {
    --rw-badge-background: var(--rw-badge-dark-background, #ffffff);
    --rw-badge-color: var(--rw-badge-dark-color, #000000);
  }

  .opacity-0 {
    opacity: 0;
  }

  .translate-x-full {
    transform: translateX(150%);
  }

  .-translate-x-full {
    transform: translateX(-150%);
  }

  .left {
    left: var(--rw-badge-position-offset, 20px);
  }

  .right {
    right: var(--rw-badge-position-offset, 20px);
  }

  .top {
    top: var(--rw-badge-position-offset, 20px);
  }

  .bottom {
    bottom: var(--rw-badge-position-offset, 20px);
  }

  #rapidweaverBadge {
    background: var(--rw-badge-background, #000000);
    border-radius: var(--rw-badge-border-radius, 6px);
    box-shadow: var(--rw-badge-box-shadow, 0 1px 5px rgba(0,0,0,0.10), 0 1px 2px rgba(0,0,0,0.2));
    padding: var(--rw-badge-padding, 6px 10px);
    position: var(--rw-badge-position, fixed);
    transition: all var(--rw-badge-transition-duration, 300ms) ease-in-out;
    z-index: var(--rw-badge-z-index, 99999);
    text-decoration: var(--rw-badge-text-decoration, none) !important;
  }

  #rapidweaverBadgeContent,
  #rapidweaverBadge a:link,
  #rapidweaverBadge a:visited,
  #rapidweaverBadge a:hover,
  #rapidweaverBadge a:active {
    text-decoration: var(--rw-badge-text-decoration, none) !important;
    align-items: var(--rw-badge-align-items, center) !important;
    display: var(--rw-badge-display, flex) !important;
    gap: var(--rw-badge-gap, 6px) !important;
  }

  ::slotted(img) {
    height: var(--rw-badge-image-height, 28px) !important;
    width: var(--rw-badge-image-width, 28px) !important;
  }

  ::slotted(*) {
    color: var(--rw-badge-color, #ffffff) !important;
    display: block !important;
    font-family: var(--rw-badge-font-family, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol") !important;
    font-size: var(--rw-badge-font-size, 12px) !important;
    font-weight: var(--rw-badge-font-weight, bold) !important;
    margin:0 !important;
    padding: 0 !important;
    text-transform: var(--rw-badge-text-transform, none) !important;
    text-decoration: var(--rw-badge-text-decoration, none) !important;
  }
</style>

<div id="rapidweaverBadge">
  <div id="rapidweaverBadgeContent">
    <slot></slot>
  </div>
</div>
`;

class RapidWeaverBadge extends HTMLElement {

  $badge;
  $badgeContent;
  $defaults = {
    positionX: 'left',
    positionY: 'bottom',
    mode: 'auto',
    transition: 'slide',
    delayType: 'time',
    delay: 1000,
    url: null,
    target: '_self',
    hideOn: null
  };

  constructor() {
    super();
    this.attachShadow({ mode: 'open' });
  }

  connectedCallback() {
    this.shadowRoot.appendChild(RapidWeaverBadgeTemplate.content.cloneNode(true));
    this.$badge = this.shadowRoot.getElementById('rapidweaverBadge')
    this.$badgeContent = this.shadowRoot.getElementById('rapidweaverBadgeContent')

    if (this.url && this.url !== '') {
      this.addLink();
    }

    this.addTransitions();
    this.setPosition();
    this.setTheme();
    this.addHiddenStyles();
    window.matchMedia("(prefers-color-scheme: dark)")
      .addEventListener("change",   e => this.setTheme());
  }

  /**
   * @returns {string}
   * @private
   * @memberof RapidWeaverBadge
   * @description Work out when to show/hide the badge
  **/
  addHiddenStyles() {
    const styleNode = Array.from(this.shadowRoot.childNodes)
      .filter(childNode => childNode.nodeName === 'STYLE');

    if (Array.isArray(this.hideWhen)) {
      this.hideWhen.forEach(hideWhen => {
        styleNode[0].innerHTML += `
          @media ${hideWhen} {
            #rapidweaverBadge {
              display: none !important;
              visibility: hidden !important;
            }
          }
        `;
      });
    }

  }

  /**
   * @returns {string}
   * @private
   * @memberof RapidWeaverBadge
   * @description Wrap the content in a link
  **/
  addLink() {
    const link = document.createElement('a');
    link.href = this.url;
    link.id = "rapidweaverBadgeContent";
    link.target = this.target;
    link.innerHTML = this.$badgeContent.innerHTML;

    this.$badge.innerHTML = '';
    this.$badge.appendChild(link);

    return link;
  }

  /**
   * @returns {null}
   * @private
   * @memberof RapidWeaverBadge
   * @description Add the transitions to the badge
  **/
  addTransitions() {
    if (this.transition === 'none' || !this.transition || !this.delayType) {
      return;
    }

    this.$badge.classList.add(this.transitionClasses);

    // if the delayType is 'time' them we add a timeout
    if (this.delayType === 'time') {
      setTimeout(() => {
        this.$badge.classList.remove(this.transitionClasses);
      }, this.delay);
    }

    // if the delayType is 'scroll' we add a scroll listener
    if (this.delayType === 'scroll') {
      window.addEventListener('scroll', () => {
        const scrollTop = window.pageYOffset || document.documentElement.scrollTop;
        const scrollPercent = (scrollTop / (document.body.scrollHeight - window.innerHeight)) * 100;

        if (scrollPercent >= this.delay) {
          this.$badge.classList.remove(this.transitionClasses);
        }
      });
    }
  }

  /**
   * @returns {string}
   * @private
   * @memberof RapidWeaverBadge
   * @description Set the theme to light, dark, or auto
  **/
  setTheme() {;
    if (this.mode === 'auto') {
      const theme = window.matchMedia("(prefers-color-scheme: dark)").matches ? 'dark' : 'light';
      this.$badge.setAttribute('data-theme', theme);
      return;
    }

    this.$badge.setAttribute('data-theme', this.mode);
  }

  /**
   * @returns {null}
   * @private
   * @memberof RapidWeaverBadge
   * @description Set the position to top or bottom, and left or right
  **/
  setPosition() {
    const classes = {
      'left': 'left',
      'right': 'right',
      'top': 'top',
      'bottom': 'bottom',
    }

    this.$badge.classList.add(classes[this.positionX]);
    this.$badge.classList.add(classes[this.positionY]);
  }

  /**
   * @returns {string}
   * @private
   * @memberof RapidWeaverBadge
   * @description Set the transition classes
  **/
  get transitionClasses() {
    const name = this.transition === 'slide' ? `slide-${this.positionX}` : this.transition;
    return {
      'none': '',
      'fade': 'opacity-0',
      'slide-left': '-translate-x-full',
      'slide-right': 'translate-x-full'
    }[name];
  }

  /**
   * @returns {string}
   * @private
   * @memberof RapidWeaverBadge
   * @description Get the URL attribute
  **/
  get url() {
    return this.getAttribute("url");
  }

  /**
   * @returns {string}
   * @private
   * @memberof RapidWeaverBadge
   * @description Get the delayType attribute
  **/
  get delayType() {
    const delayType = this.getAttribute("delay-type");

    if (['time', 'scroll', 'none'].indexOf(delayType) === -1) {
      return this.$defaults.delayType;
    }

    if (delayType === 'none') {
      return false;
    }

    return !delayType ? this.$defaults.delayType : delayType;
  }

  /**
   * @returns {string}
   * @private
   * @memberof RapidWeaverBadge
   * @description Get the delay attribute
  **/
  get delay() {
    const delay = parseInt(this.getAttribute("delay"));
    return !delay ? this.$defaults.delay : delay;
  }

  /**
   * @returns {string}
   * @private
   * @memberof RapidWeaverBadge
   * @description Get the transition attribute
  **/
  get transition() {
    const transition = this.getAttribute("transition");

    if (['none', 'fade', 'slide'].indexOf(transition) === -1) {
      return this.$defaults.transition;
    }

    if (transition === 'none') {
      return false;
    }

    return !transition ? this.$defaults.transition : transition;
  }

  /**
   * @returns {string}
   * @private
   * @memberof RapidWeaverBadge
   * @description Get the position-x attribute
  **/
  get positionX() {
    const positionX = this.getAttribute("position-x");

    if (['left', 'right'].indexOf(positionX) === -1) {
      return this.$defaults.positionX;
    }

    return !positionX ? this.$defaults.positionX : positionX;
  }

  /**
   * @returns {string}
   * @private
   * @memberof RapidWeaverBadge
   * @description Get the position-y attribute
  **/
  get positionY() {
    const positionY = this.getAttribute("position-y");

    if (['top', 'bottom'].indexOf(positionY) === -1) {
      return this.$defaults.positionY;
    }

    return !positionY ? this.$defaults.positionY : positionY;
  }

  /**
   * @returns {string}
   * @private
   * @memberof RapidWeaverBadge
   * @description Get the mode attribute
  **/
  get mode() {
    const mode = this.getAttribute("mode");

    if (['auto', 'light', 'dark'].indexOf(mode) === -1) {
      return this.$defaults.mode;
    }

    return !mode ? this.$defaults.mode : mode;
  }

  /**
   * @returns {string}
   * @private
   * @memberof RapidWeaverBadge
   * @description Get the target attribute
  **/
  get target() {
    const target = this.getAttribute("target");
    return !target ? this.$defaults.target : target;
  }

  /**
   * @returns {array}
   * @private
   * @memberof RapidWeaverBadge
   * @description Get the hide-on attribute
  **/
   get hideWhen() {
    const hideWhen = this.getAttribute("hide-when");
    return !hideWhen ? this.$defaults.hideWhen : hideWhen.split(',');
  }
}

window.customElements.define('rapidweaver-badge', RapidWeaverBadge);
